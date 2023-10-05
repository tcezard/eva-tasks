package eva3399

import eva3399.ClashingSSHashes
import org.slf4j.LoggerFactory
import org.springframework.data.mongodb.core.query.Update
import uk.ac.ebi.ampt2d.commons.accession.core.models.EventType
import uk.ac.ebi.eva.accession.core.EVAObjectModelUtils
import uk.ac.ebi.eva.accession.core.model.dbsnp.DbsnpSubmittedVariantOperationEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantInactiveEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantOperationEntity
import uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment
import uk.ac.ebi.eva.groovy.commons.RetryableBatchingCursor
import uk.ac.ebi.eva.remapping.ingest.batch.io.SubmittedVariantDiscardPolicy

import static org.springframework.data.mongodb.core.query.Criteria.where
import static org.springframework.data.mongodb.core.query.Query.query
import static uk.ac.ebi.eva.groovy.commons.EVADatabaseEnvironment.*
import static uk.ac.ebi.eva.remapping.ingest.batch.io.SubmittedVariantDiscardPolicy.SubmittedVariantDiscardDeterminants

class RemediateIndels {
    static def logger = LoggerFactory.getLogger(RemediateIndels.class)
    String assembly
    EVADatabaseEnvironment prodEnv
    EVADatabaseEnvironment devEnv
    boolean isTargetAssembly

    RemediateIndels() {}

    RemediateIndels(String assembly, EVADatabaseEnvironment prodEnv, EVADatabaseEnvironment devEnv,
                    boolean isTargetAssembly=false) {
        this.assembly = assembly
        this.prodEnv = prodEnv
        this.devEnv = devEnv
        this.isTargetAssembly = isTargetAssembly
    }

    // Remove old SS hash tacked on the start of remapping ID. While this will be useful later on, for now we don't need it!
    // See here: https://github.com/EBIvariation/eva-tasks/blob/27e4a1bda59fb30e8baabd2500119ae9ffacf05d/tasks/eva_3371/src/main/groovy/eva3371/eva3371_detect_unnormalized_indels.groovy#L169
    def _undoEVA3371Hack = { SubmittedVariantEntity sve ->
        if (Objects.isNull(sve.remappedFrom)) sve.setRemappingId(null)
        if (Objects.nonNull(sve.remappingId)) {
            if (sve.remappingId.contains("_")) {
                sve.setRemappingId(sve.remappingId.split("_")[1])
            }
        }
        return sve
    }

    // Get SS which are candidates for merging - collected in EVA-3380
    def _getMergeableSS = {List<SubmittedVariantEntity> sves ->
        def ssHashesToFind = sves.collect{it.hashedMessage}
        def clashingHashRecords = this.devEnv.mongoTemplate.find(query(where("_id").in(ssHashesToFind)),
                ClashingSSHashes.class)
        def allAccessionsInClashingRecords = clashingHashRecords
                .collect{it.clashingSS.collect{it.accession}}.flatten()
        def existingAccessionsInProd = [sveClass, dbsnpSveClass].collect{collectionClass ->
            this.prodEnv.mongoTemplate.find(query(where("accession").in(allAccessionsInClashingRecords)
                    .and("seq").is(this.assembly)), collectionClass)
        }.flatten().collect{it.accession}.toSet()
        clashingHashRecords.each{clashingHashRecord ->
            // Invalidate SS which no longer exist due to deletions elsewhere in this script
            clashingHashRecord.clashingSS.removeIf{!(existingAccessionsInProd.contains(it.accession))}
            clashingHashRecord.clashingSS = clashingHashRecord.clashingSS.collect{_undoEVA3371Hack(it)}
        }
        return clashingHashRecords
    }

    def _prioritise = {SubmittedVariantEntity sve1, SubmittedVariantEntity sve2 ->
        def sve1DD = new SubmittedVariantDiscardDeterminants(sve1, sve1.accession, sve1.remappedFrom, sve1.createdDate)
        def sve2DD = new SubmittedVariantDiscardDeterminants(sve2, sve2.accession, sve2.remappedFrom, sve2.createdDate)
        return SubmittedVariantDiscardPolicy.prioritise(sve1DD, sve2DD).sveToKeep.sve
    }

    // Given a set of merge candidates in a ClashingSSHashes object, determine the target SS
    def _getMergeTargetAndMergees = { ClashingSSHashes clashingSSHashRecord ->
        // If additional merge targets were deleted elsewhere in this script for some reason
        // just return the merge target with no mergees
        if (clashingSSHashRecord.clashingSS.size() == 1) return [clashingSSHashRecord.clashingSS[0], []]
        def mergeTarget = clashingSSHashRecord.clashingSS[0]
        (1..(clashingSSHashRecord.clashingSS.size()-1)).each{i ->
            mergeTarget = _prioritise(mergeTarget, clashingSSHashRecord.clashingSS[i])
        }
        return [mergeTarget,
                clashingSSHashRecord.clashingSS.findAll{it.accession != mergeTarget.accession}]
    }

    def _getSSRecordsBeforeNorm = {List<SubmittedVariantEntity> svesToNormalize ->
        return [sveClass, dbsnpSveClass].collect {collectionClass ->
            this.prodEnv.mongoTemplate.find(query(where("accession").in(svesToNormalize.collect{it.accession})
                    .and("seq").is(this.assembly)), collectionClass)
        }.flatten()
    }

    def _writeSSLocusUpdateOps = {List<SubmittedVariantEntity> svesWithUpdatedLocus,
                                  List<SubmittedVariantEntity> ssRecordsBeforeNorm ->
        def (evaSVOEOps, dbsnpSVOEOps) = [new ArrayList<>(), new ArrayList<>()]
        def oldSSRecordsGroupedBySSID = ssRecordsBeforeNorm.groupBy { it.accession }
        svesWithUpdatedLocus.each {sve ->
            if (oldSSRecordsGroupedBySSID.containsKey(sve.accession)) {
                def ssUpdatedOp = new SubmittedVariantOperationEntity()
                ssUpdatedOp.fill(EventType.UPDATED, sve.accession, "EVA3399 - SS received new locus after normalization",
                        oldSSRecordsGroupedBySSID.get(sve.accession).collect { new SubmittedVariantInactiveEntity(it) })
                ssUpdatedOp.setId("EVA3399_UPD_LOCUS_${this.assembly}_${sve.accession}_${sve.hashedMessage}".toString())
                (sve.accession >= 5e9) ? evaSVOEOps.add(ssUpdatedOp) : dbsnpSVOEOps.add(ssUpdatedOp)
            }
        }
        this.prodEnv.bulkInsertIgnoreDuplicates(evaSVOEOps, svoeClass)
        this.prodEnv.bulkInsertIgnoreDuplicates(dbsnpSVOEOps, dbsnpSvoeClass)
    }

    def _writeMergeTarget = {SubmittedVariantEntity mergeTarget ->
        def singletonListWithMergeTarget = Collections.singletonList(mergeTarget)
        def mergeTargetCollection = (mergeTarget.accession >= 5e9) ? sveClass : dbsnpSveClass
        def hashExistsInProd = ([sveClass, dbsnpSveClass].collect {collectionClass ->
            this.prodEnv.mongoTemplate.find(query(where("_id").is(mergeTarget.hashedMessage)), collectionClass)
        }.flatten().size() > 0)
        def mergeTargetRecordsBeforeNorm = _getSSRecordsBeforeNorm(singletonListWithMergeTarget)
        // The following replaces the existing SS entry with the hash with the mergeTarget record
        // Or creates a new entry if no current entry exists
        this.prodEnv.mongoTemplate.save(mergeTarget, prodEnv.mongoTemplate.getCollectionName(mergeTargetCollection))
        // If the normalized hash did not previously exist in PROD, it means we have updated the locus of an existing SS
        // Write an operation to indicate that the locus was updated for the SS
        if (!hashExistsInProd) _writeSSLocusUpdateOps(singletonListWithMergeTarget, mergeTargetRecordsBeforeNorm)
    }

    def _removeSves = {List<SubmittedVariantEntity> sves, boolean deleteRemapped=false ->
        def oldSSIDs = sves.collect{it.accession}
        [sveClass, dbsnpSveClass].each{collectionClass ->
            // Remove all occurrences of this SS including in remapped assemblies
            this.prodEnv.mongoTemplate.findAllAndRemove(query(where("accession").in(oldSSIDs)
                    .and("seq").is(this.assembly)), collectionClass)
            if (deleteRemapped) {
                this.prodEnv.mongoTemplate.find(query(where("accession").in(oldSSIDs)
                        .and("remappedFrom").is(this.assembly)), collectionClass).each {sve ->
                    def opClass = collectionClass.equals(sveClass)?
                            new SubmittedVariantOperationEntity(): new DbsnpSubmittedVariantOperationEntity()
                    opClass.fill(EventType.DISCARDED, "EVA3399 - SS was merged in the source assembly",
                            Arrays.asList(new SubmittedVariantInactiveEntity(sve)))
                    opClass.setId("EVA3399_DISCARDED_" +
                            "${sve.referenceSequenceAccession}_${sve.accession}_${sve.hashedMessage}".toString())
                    this.prodEnv.mongoTemplate.save(opClass)
                }
                this.prodEnv.mongoTemplate.findAllAndRemove(query(where("accession").in(oldSSIDs)
                        .and("remappedFrom").is(this.assembly)), collectionClass)
            }
        }
    }

    def _writeMergeOps = {List<SubmittedVariantEntity> mergees, SubmittedVariantEntity mergeTarget ->
        mergees.each { mergee ->
            logger.info("Merging ${mergee.accession} to ${mergeTarget.accession}...")
            def mergeOpId = "EVA3399_MERGED_${this.assembly}_${mergee.accession}_${mergeTarget.accession}".toString()
            def mergeOpClass = (mergee.accession >= 5e9) ? svoeClass : dbsnpSvoeClass
            def mergeOp = new SubmittedVariantOperationEntity()
            mergeOp.fill(EventType.MERGED, mergee.accession, mergeTarget.accession,
                    "EVA3399 - Hash collision with another SS after normalization",
                    Arrays.asList(new SubmittedVariantInactiveEntity(mergee)))
            mergeOp.setId(mergeOpId)
            this.prodEnv.mongoTemplate.save(mergeOp, this.prodEnv.mongoTemplate.getCollectionName(mergeOpClass))
        }
    }

    def _updateSSWithRS = {List<SubmittedVariantEntity> svesToUpdate, Map<String, Long> rsIDsByHash ->
        _clearExistingRs(svesToUpdate)
        def svesGroupedByRsHash = svesToUpdate.groupBy { EVAObjectModelUtils.getClusteredVariantHash(it)}
        rsIDsByHash.each {rsHash, rsID ->
            def sves = svesGroupedByRsHash.getOrDefault(rsHash, [])
            def svesHash = sves.collect{it.hashedMessage}
            sves.each {it.setClusteredVariantAccession(rsID)}
            sves.each{logger.info("Updating SS ${it.accession} with RS ${rsID}...")}
            [sveClass, dbsnpSveClass].each {collectionClass ->
                this.prodEnv.mongoTemplate.updateMulti(query(where("_id").in(svesHash)),
                        Update.update("rs", rsID), collectionClass)
            }
        }
    }

    def _clearExistingRs = {List<SubmittedVariantEntity> sves ->
        def svesHash = sves.collect{it.hashedMessage}
        sves.each{it.setClusteredVariantAccession(null)}
        [sveClass, dbsnpSveClass].each {ssClass ->
            this.prodEnv.mongoTemplate.updateMulti(query(where("_id").in(svesHash)), new Update().unset("rs"), ssClass)
        }
    }

    def _assignNewRs = {List<SubmittedVariantEntity> svesToUpdate ->
        def rsRecordsToCreate = svesToUpdate.collect { EVAObjectModelUtils.toClusteredVariant(it)}
        def rsIDsByHash = this.prodEnv.clusteredVariantAccessioningService.getOrCreate(rsRecordsToCreate).collectEntries { [it.hash, it.accession] }
        _updateSSWithRS(svesToUpdate, rsIDsByHash)
    }

    def _assignRs = {List<SubmittedVariantEntity> svesToUpdate ->
        def svesToBeAssignedExistingRS = svesToUpdate.findAll{Objects.nonNull(it.clusteredVariantAccession)}
        def existingRsHashesToFind = svesToUpdate.collect{EVAObjectModelUtils.getClusteredVariantHash(it)}.toSet()
        def existingRsIDsByHash = [cveClass, dbsnpCveClass].collect{collectionClass ->
            this.prodEnv.mongoTemplate.find(query(where("_id")
                    .in(existingRsHashesToFind)), collectionClass)
        }.flatten().collectEntries{[it.hashedMessage, it.accession]}
        _updateSSWithRS(svesToBeAssignedExistingRS, existingRsIDsByHash)
        if (this.isTargetAssembly) {
            _assignNewRs(svesToUpdate.findAll {
                !(existingRsIDsByHash.containsKey(EVAObjectModelUtils.getClusteredVariantHash(it)))})
        }
    }

    def _merge = {SubmittedVariantEntity mergeTarget, List<SubmittedVariantEntity> mergees ->
        _writeMergeTarget(mergeTarget)
        // Remove mergees in this and remapped assemblies as well
        _removeSves(mergees, true)
        _assignRs(Collections.singletonList(mergeTarget))
        _writeMergeOps(mergees, mergeTarget)
    }

    def mergeInProdEnv = { List<ClashingSSHashes> clashingHashes ->
        def mergeTargetsAndMergeeList = clashingHashes.collect { _getMergeTargetAndMergees(it) }
        mergeTargetsAndMergeeList.each {mergeTargetAndMergees ->
            def (mergeTarget, mergees) = mergeTargetAndMergees
            _merge(mergeTarget, mergees)
        }
    }

    def insertIntoProdEnv = {List<SubmittedVariantEntity> svesToInsert ->
        def correspondingOldSveRecords = _getSSRecordsBeforeNorm(svesToInsert)
        _removeSves(svesToInsert)
        this.prodEnv.mongoTemplate.insert(svesToInsert.findAll{it.accession >= 5e9}, sveClass)
        this.prodEnv.mongoTemplate.insert(svesToInsert.findAll{it.accession < 5e9}, dbsnpSveClass)
        _writeSSLocusUpdateOps(svesToInsert, correspondingOldSveRecords)
        _assignRs(svesToInsert)
    }

    def runRemediation = {
        def evaAndDbsnpSveCursorsDev = [sveClass, dbsnpSveClass].collect { collectionClass ->
            new RetryableBatchingCursor<>(where("seq").is(assembly), devEnv.mongoTemplate, collectionClass)
        }
        def numRecordsProcessed = 0
        evaAndDbsnpSveCursorsDev.each {
            it.each {List<SubmittedVariantEntity> sves ->
                sves = sves.collect{_undoEVA3371Hack(it)}
                def mergeableSS = _getMergeableSS(sves)
                def mergeableSSHashes = mergeableSS.collect{it.ssHash}.toSet()
                def unmergeableSS = sves.findAll{!(mergeableSSHashes.contains(it.hashedMessage))}
                mergeInProdEnv(mergeableSS)
                insertIntoProdEnv(unmergeableSS)
                numRecordsProcessed += sves.size()
                logger.info("${numRecordsProcessed} SS processed so far...")
            }
        }
    }
}
