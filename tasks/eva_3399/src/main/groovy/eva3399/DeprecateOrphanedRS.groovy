package eva3399

import uk.ac.ebi.ampt2d.commons.accession.core.models.EventType
import uk.ac.ebi.eva.accession.core.EVAObjectModelUtils
import uk.ac.ebi.eva.accession.core.batch.io.ClusteredVariantDeprecationWriter
import uk.ac.ebi.eva.accession.core.model.dbsnp.DbsnpClusteredVariantOperationEntity
import uk.ac.ebi.eva.accession.core.model.eva.ClusteredVariantOperationEntity
import uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment
import uk.ac.ebi.eva.groovy.commons.RetryableBatchingCursor

import static org.springframework.data.mongodb.core.query.Criteria.where
import static org.springframework.data.mongodb.core.query.Query.query
import static uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment.*

class DeprecateOrphanedRS {
    String assembly
    EVADatabaseEnvironment dbEnv

    DeprecateOrphanedRS(String assembly, EVADatabaseEnvironment dbEnv) {
        this.assembly = assembly
        this.dbEnv = dbEnv
    }

    void deprecate() {
        def deprecationIdSuffix = "EVA3399"
        def rsDeprecationWriter = new ClusteredVariantDeprecationWriter(this.assembly, this.dbEnv.mongoTemplate,
                this.dbEnv.submittedVariantAccessioningService,
                this.dbEnv.springApplicationContext.getBean("accessioningMonotonicInitRs", Long.class),
                deprecationIdSuffix, "RS orphaned due to locus not used by any SS")
        def opsOfImpactedSS = [svoeClass, dbsnpSvoeClass].collect{new RetryableBatchingCursor<>(
                where("_id").regex("^EVA3399_UPD_LOCUS_${this.assembly}_.*"), this.dbEnv.mongoTemplate, it)}
        opsOfImpactedSS.each{it.each{opsInBatch ->
            def ssOpsGroupedByRSHashes = opsInBatch.groupBy{
                EVAObjectModelUtils.getClusteredVariantHash(it.inactiveObjects[0])}
            def rsToDeprecate = [cveClass, dbsnpCveClass].collect{
                this.dbEnv.mongoTemplate.find(query(where("_id").in(ssOpsGroupedByRSHashes.keySet())), it)
            }.flatten()
            rsDeprecationWriter.write(rsToDeprecate)
            def deprecationOpIds = ssOpsGroupedByRSHashes.keySet().collect{
                "RS_DEPRECATED_${deprecationIdSuffix}_${it}".toString()}
            // Write RS merge operation record for RS that ended up getting deprecated above
            [cvoeClass, dbsnpCvoeClass].each {opClass ->
                this.dbEnv.mongoTemplate.find(query(where("_id").in(deprecationOpIds)), opClass).each {op ->
                    def newSSHashToLookup =
                            ssOpsGroupedByRSHashes[op.inactiveObjects[0].hashedMessage][0].id.split("_")[-1]
                    def oldRS =
                            ssOpsGroupedByRSHashes[op.inactiveObjects[0].hashedMessage][0].inactiveObjects[0].clusteredVariantAccession
                    def newRSForThisSS = [sveClass, dbsnpSveClass].collect{ssClass ->
                        this.dbEnv.mongoTemplate.find(query(where("_id").is(newSSHashToLookup)
                                .and("rs").ne(oldRS)), ssClass)
                    }.flatten().collect{it.clusteredVariantAccession}
                    if (newRSForThisSS.size() > 0) {
                        def mergeOpToWrite = opClass.equals(cvoeClass) ?
                                new ClusteredVariantOperationEntity(): new DbsnpClusteredVariantOperationEntity()
                        mergeOpToWrite.fill(EventType.MERGED, op.accession, newRSForThisSS[0],
                                "EVA3399 - SS with updated locus received new RS", op.inactiveObjects)
                        mergeOpToWrite.setId("EVA3399_MERGED_${this.assembly}_${op.accession}_${newRSForThisSS[0]}")
                        this.dbEnv.mongoTemplate.save(mergeOpToWrite)
                    }
                }
            }
            // We also need to remove the deprecation ops so that they don't show up in the website
            [cvoeClass, dbsnpCvoeClass].each {opClass ->
                this.dbEnv.mongoTemplate.findAllAndRemove(query(where("_id").in(deprecationOpIds)), opClass)}
        }}
    }
}
