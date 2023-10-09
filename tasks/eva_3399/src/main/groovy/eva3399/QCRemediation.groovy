package eva3399

import org.slf4j.Logger
import org.slf4j.LoggerFactory
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantOperationEntity
import uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment
import uk.ac.ebi.eva.groovy.commons.RetryableBatchingCursor

import static org.springframework.data.mongodb.core.query.Criteria.where
import static org.springframework.data.mongodb.core.query.Query.query
import static uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment.*

class QCRemediation {
    static Logger logger = LoggerFactory.getLogger(QCRemediation.class)
    String assemblyToQC
    EVADatabaseEnvironment prodEnv
    EVADatabaseEnvironment devEnv
    boolean isTargetAssembly
    private RemediateIndels remediateIndelsObj

    QCRemediation(String assemblyToQC, EVADatabaseEnvironment prodEnv, EVADatabaseEnvironment devEnv,
                  boolean isTargetAssembly) {
        this.assemblyToQC = assemblyToQC
        this.prodEnv = prodEnv
        this.devEnv = devEnv
        this.isTargetAssembly = isTargetAssembly
        this.remediateIndelsObj = new RemediateIndels(this.assemblyToQC, this.prodEnv, this.devEnv,
                this.isTargetAssembly)
    }

    void checkMergedSvesInProd() {
        def ssOpCursors = [svoeClass, dbsnpSvoeClass].collect{opClass ->
            new RetryableBatchingCursor<>(where("_id").regex("EVA3399_MERGED_${this.assemblyToQC}_.*"),
                    this.prodEnv.mongoTemplate, opClass)
        }
        ssOpCursors.each {it.each{List<SubmittedVariantOperationEntity> svoes ->
            def expectedMergeeSveIDs = svoes.collect { it.inactiveObjects[0].accession }.toSet()
            // Ensure that mergee SVEs are not present in PROD
            ensureMergeeSvesNotPresent(expectedMergeeSveIDs)
            // Ensure that merged SVEs obey prioritization
            ensureCorrectMergeTarget(svoes)
        }}
    }

    void ensureProperRSAssignment() {
        def ssOpCursors = [svoeClass, dbsnpSvoeClass].collect{opClass ->
            new RetryableBatchingCursor<>(where("_id").regex("EVA3399_UPD_LOCUS_${this.assemblyToQC}_.*"),
                    this.prodEnv.mongoTemplate, opClass)
        }
        ssOpCursors.each {it.each { List<SubmittedVariantOperationEntity> svoes ->
            def ssHashesToLookup = svoes.findAll {
                Objects.isNull(it.inactiveObjects[0].clusteredVariantAccession)
            }.collect { it.inactiveObjects.collect { it.hashedMessage } }.flatten()
            [sveClass, dbsnpSveClass].each { ssClass ->
                this.prodEnv.mongoTemplate.find(query(where("_id").in(ssHashesToLookup).and("rs").exists(true)),
                        ssClass).each {
                    logger.error("SS with hash ${it.hashedMessage} in ${this.prodEnv.mongoTemplate.getCollectionName(ssClass)} " +
                            "should not have been assigned an RS because it did not have an RS pre-remediation!")
                }
            }
        }}
    }

    void ensureCorrectMergeTarget(List<SubmittedVariantOperationEntity> svoes) {
        def svoesGroupedByHash = svoes.groupBy{it.inactiveObjects[0].hashedMessage}
        def mergeTargetsInProdGroupedByHash = [sveClass, dbsnpSveClass].collect{
            this.prodEnv.mongoTemplate.find(query(where("_id").in(svoesGroupedByHash.keySet())), it)
        }.flatten().collectEntries {[it.hashedMessage, it]}
        (svoesGroupedByHash.keySet() - mergeTargetsInProdGroupedByHash.keySet()).each {
            logger.error("SS hash ${it} not found in PROD!")
        }
        mergeTargetsInProdGroupedByHash.each {ssHash, prodSve ->
            def mergeeSves = svoesGroupedByHash.get(ssHash).collect{it.inactiveObjects}.flatten().collect{it.toSubmittedVariantEntity()}
            mergeeSves.each {mergeeSve ->
                def expectedID = this.remediateIndelsObj._prioritise(mergeeSve, prodSve).accession
                if (expectedID != prodSve.accession) {
                    logger.error("Submitted variant with ${prodSve.accession} in PROD " +
                            "should not have been the merge target! Expected ${expectedID}!")
                }
            }
        }
    }

    void ensureMergeeSvesNotPresent(Set<Long> expectedMergeeSveIDs) {
        [sveClass, dbsnpSveClass].each { ssClass ->
            this.prodEnv.mongoTemplate.find(query(where("accession").in(expectedMergeeSveIDs)
                    .and("seq").is(this.assemblyToQC)), ssClass).each { SubmittedVariantEntity sve ->
                logger.error("Submitted variant with SS ID ${sve.accession} should not be " +
                        "present in ${this.prodEnv.mongoTemplate.getCollectionName(ssClass)}!")
            }
        }
    }

    void prodQCChecks(List<SubmittedVariantEntity> devSves) {
        def devSveHashes = devSves.collect{it.hashedMessage}
        def prodSves = [sveClass, dbsnpSveClass].collect{
            this.prodEnv.mongoTemplate.find(query(where("_id").in(devSveHashes)), it)
        }.flatten()
        checkEquivalentSves(devSves, prodSves)
        if (this.isTargetAssembly) checkBackPropRS(prodSves)
    }

    static boolean _areSvesEqual(SubmittedVariantEntity sve1, SubmittedVariantEntity sve2) {
        return sve1.taxonomyAccession == sve2.taxonomyAccession && sve1.start == sve2.start &&
                Objects.equals(sve1.referenceSequenceAccession, sve2.referenceSequenceAccession) &&
                Objects.equals(sve1.projectAccession, sve2.projectAccession) &&
                Objects.equals(sve1.contig, sve2.contig) &&
                Objects.equals(sve1.referenceAllele, sve2.referenceAllele) &&
                Objects.equals(sve1.alternateAllele, sve2.alternateAllele) &&
                Objects.equals(sve1.isSupportedByEvidence(), sve2.isSupportedByEvidence()) &&
                Objects.equals(sve1.isAssemblyMatch(), sve2.isAssemblyMatch()) &&
                Objects.equals(sve1.isAllelesMatch(), sve2.isAllelesMatch()) &&
                Objects.equals(sve1.isValidated(), sve2.isValidated()) &&
                Objects.equals(sve1.remappedFrom, sve2.remappedFrom) &&
                Objects.equals(sve1.remappedDate, sve2.remappedDate) &&
                Objects.equals(sve1.mapWeight, sve2.mapWeight)
    }

    static void checkEquivalentSves(List<SubmittedVariantEntity> devSves, List<SubmittedVariantEntity> prodSves) {
        def (devSvesGroupedByHashes, prodSvesGroupedByHashes) =
        [devSves, prodSves].collect { it.collectEntries {[it.hashedMessage, it]} }
        devSvesGroupedByHashes.each {ssHash, devSve ->
            if (!prodSvesGroupedByHashes.containsKey(ssHash)) {
                logger.error("Submitted Variant with hash ${ssHash} not found in PROD!")
            } else {
                def prodSve = prodSvesGroupedByHashes[ssHash]
                if (prodSve.accession != devSve.accession) {
                    logger.error("Submitted Variant with hash ${ssHash} is not associated " +
                            "with SS ID ${devSve.accession} in PROD!")
                }

                if (!_areSvesEqual(devSve, prodSve)) {
                    logger.error("Submitted Variant in PROD ${prodSve} is not equivalent to " +
                            "the expected submitted variant ${devSve}!")
                }
            }
        }
    }

    void checkBackPropRS(List<SubmittedVariantEntity> prodSves) {
        def remappedSves = prodSves.findAll{
            Objects.nonNull(it.remappedFrom) && Objects.nonNull(it.clusteredVariantAccession)}
        Map<Long, SubmittedVariantEntity> remappedSvesGroupedByAccession =
                remappedSves.collectEntries {[it.accession, it]}
        def sourceAssemblies = remappedSves.collect{it.remappedFrom}.toSet()
        def remappedSveIDs = prodSves.collect{it.accession}
        def remappedSveRSIDs = prodSves.collect{it.clusteredVariantAccession}.toSet()
        List<SubmittedVariantEntity> svesWithBackPropRSInProd = [sveClass, dbsnpSveClass].collect {
            this.prodEnv.mongoTemplate.find(query(where("accession").in(remappedSveIDs)
                    .and("seq").in(sourceAssemblies).and("backPropRS").in(remappedSveRSIDs)), it)}.flatten()
        (remappedSveRSIDs - svesWithBackPropRSInProd.collect{it.backPropagatedVariantAccession}.toSet()).each {
            logger.error("No backpropagated RS entry for RS ${it} in the source assembly " +
                    "originating from the remapped assembly ${this.assemblyToQC}!")
        }
        svesWithBackPropRSInProd.each {sve ->
            def expectedBackPropRS = remappedSvesGroupedByAccession[sve.accession].clusteredVariantAccession
            if (sve.backPropagatedVariantAccession != expectedBackPropRS) {
                logger.error("Expected ${expectedBackPropRS} " +
                        "but found ${sve.backPropagatedVariantAccession} for SS ${sve}!")
            }
        }
    }

    void runQC() {
        // Ensure that all unmergeable SS hashes are assigned to the proper accessions
        [sveClass, dbsnpSveClass].each{ssClass ->
            new RetryableBatchingCursor<>(where("seq").is(this.assemblyToQC), this.devEnv.mongoTemplate, ssClass).each{sves ->
                sves = sves.collect{this.remediateIndelsObj._undoEVA3371Hack(it)}
                def ssHashesToFind = sves.collect{it.hashedMessage}
                List<ClashingSSHashes> mergeableSS = this.devEnv.mongoTemplate.find(query(where("_id")
                        .in(ssHashesToFind)), ClashingSSHashes.class)
                def mergeableSSHashes = mergeableSS.collect{it.ssHash}.toSet()
                def unmergeableSS = sves.findAll{!(mergeableSSHashes.contains(it.hashedMessage))}
                prodQCChecks(unmergeableSS)

                def mergeTargetAndMergeeList = mergeableSS.collect {
                    this.remediateIndelsObj._getMergeTargetAndMergees(it) }
                def mergeTargetList = mergeTargetAndMergeeList.collect{it[0]}
                def mergeeList = mergeTargetAndMergeeList.collect{it[1]}
                prodQCChecks(mergeTargetList)
            }
        }
        checkMergedSvesInProd()
        // If this is not a target assembly, ensure no existing RS is assigned for SS that did not previously have it
        if(!this.isTargetAssembly) ensureProperRSAssignment()
    }
}
