package src

@Grab(group = 'uk.ac.ebi.eva', module = 'eva-accession-clustering', version = '0.6.10-SNAPSHOT')
@Grab(group = 'uk.ac.ebi.eva', module = 'eva-accession-core', version = '0.6.10-SNAPSHOT')
@Grab(group = 'uk.ac.ebi.eva', module = 'variation-commons-batch', version = '0.8.1')

import com.mongodb.MongoCursorNotFoundException
import com.mongodb.WriteConcern
import groovy.transform.EqualsAndHashCode
import htsjdk.samtools.util.SequenceUtil
import org.apache.commons.lang3.ObjectUtils
import org.apache.commons.lang3.tuple.ImmutablePair
import org.bson.Document
import org.slf4j.Logger
import org.slf4j.LoggerFactory
import org.springframework.beans.factory.annotation.Autowired
import org.springframework.beans.factory.annotation.Value
import org.springframework.boot.CommandLineRunner
import org.springframework.boot.SpringApplication
import org.springframework.boot.autoconfigure.SpringBootApplication
import org.springframework.context.ConfigurableApplicationContext
import org.springframework.context.annotation.Import
import org.springframework.data.domain.Sort
import org.springframework.data.mongodb.core.MongoTemplate
import org.springframework.data.mongodb.core.query.Criteria
import org.springframework.data.mongodb.core.query.Query
import org.springframework.retry.annotation.Backoff
import org.springframework.retry.annotation.Retryable
import org.springframework.stereotype.Component
import uk.ac.ebi.eva.accession.clustering.configuration.InputParametersConfiguration
import uk.ac.ebi.eva.accession.clustering.configuration.batch.io.SSSplitWriterConfiguration
import uk.ac.ebi.eva.accession.clustering.configuration.batch.listeners.ListenersConfiguration
import uk.ac.ebi.eva.accession.core.configuration.nonhuman.SubmittedVariantAccessioningConfiguration
import uk.ac.ebi.eva.accession.core.model.dbsnp.DbsnpSubmittedVariantEntity
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantEntity
import uk.ac.ebi.eva.accession.core.service.nonhuman.SubmittedVariantAccessioningService
import uk.ac.ebi.eva.commons.core.models.VariantCoreFields
import uk.ac.ebi.eva.commons.core.models.VariantType
import uk.ac.ebi.eva.commons.core.models.pipeline.Variant
import uk.ac.ebi.eva.metrics.configuration.MetricConfiguration

import static org.springframework.data.mongodb.core.query.Query.query
import static org.springframework.data.mongodb.core.query.Criteria.where


@Component
@Import(value=[SSSplitWriterConfiguration.class,
        InputParametersConfiguration.class, MetricConfiguration.class, ListenersConfiguration.class,
        SubmittedVariantAccessioningConfiguration.class])
class RemediateDuplicateRemappedSS implements CommandLineRunner {

    private static final Logger logger = LoggerFactory.getLogger(RemediateDuplicateRemappedSS.class)

    @Autowired
    private SubmittedVariantAccessioningService submittedVariantAccessioningService

    @Autowired
    private MongoTemplate destinationMongoTemplate

    private MongoTemplate devMongoTemplate

    @Value('${devenv.properties}')
    private String devEnvPropertiesFile

    @Value('${assembly.to.remediate}')
    private String assemblyToRemediate

    @Autowired
    private Long accessioningMonotonicInitSs

    private Criteria criteriaToFindAssembly

    private Criteria criteriaToFindRemappedEntries

    private String remappedAttributesCollection = "remappedAttributes"

    private SSDeprecationWriter ssDeprecationWriter

    void run(String... args) {
        this.criteriaToFindAssembly = where("seq").is(assemblyToRemediate)
        this.criteriaToFindRemappedEntries = where("remappedFrom").exists(true)
        this.devMongoTemplate = EVADatabaseEnvironment.parseFrom(devEnvPropertiesFile).mongoTemplate
        this.devMongoTemplate.setWriteConcern(WriteConcern.MAJORITY)
        this.ssDeprecationWriter = new SSDeprecationWriter(this.destinationMongoTemplate, this.accessioningMonotonicInitSs)

        remediateDuplicateRemappedSSCollection(DbsnpSubmittedVariantEntity.class)
        remediateDuplicateRemappedSSCollection(SubmittedVariantEntity.class)
    }

    void remediateDuplicateRemappedSSCollection(Class<? extends SubmittedVariantEntity> collectionClassToUse) {
        int numSVEScanned = 0
        int batchIndex = 0
        String lastSeenID = null
        while(true) {
            ImmutablePair<List<? extends SubmittedVariantEntity>, String> svesAndLastSeenID =
                    getNextBatchOfDuplicateRemappedSVEs(lastSeenID, collectionClassToUse)
            if (svesAndLastSeenID != null) {
                remediateDuplicateRemappedSSNotGeneratedByRemapping(svesAndLastSeenID.left, collectionClassToUse)
                numSVEScanned += svesAndLastSeenID.left.size()
                lastSeenID = svesAndLastSeenID.right
                logger.info("Processed " + numSVEScanned + " remapped IDs in ${collectionClassToUse.getTypeName()}...")
                batchIndex += 1
            } else {
                break
            }
        }
    }

    @Retryable(value = MongoCursorNotFoundException.class, maxAttempts = 5, backoff = @Backoff(delay = 100L))
    ImmutablePair<List<? extends SubmittedVariantEntity>, String> getNextBatchOfDuplicateRemappedSVEs(
            String lastSeenID, Class<? extends SubmittedVariantEntity> collectionClassToUse) {
        String sortAttribute = "_id"
        Query queryToGetNextBatchOfDuplicateRemappedSVEs = query(this.criteriaToFindAssembly)
                .addCriteria(this.criteriaToFindRemappedEntries)
        if (Objects.nonNull(lastSeenID)) {
            queryToGetNextBatchOfDuplicateRemappedSVEs.addCriteria(where(sortAttribute).gt(lastSeenID))
        }
        queryToGetNextBatchOfDuplicateRemappedSVEs =
                queryToGetNextBatchOfDuplicateRemappedSVEs.with(Sort.by(Sort.Direction.ASC, sortAttribute)).limit(1000)
        List<? extends SubmittedVariantEntity> result =
                this.devMongoTemplate.find(queryToGetNextBatchOfDuplicateRemappedSVEs, collectionClassToUse)
        if (result.size() > 0) {
            return new ImmutablePair<>(result, result.get(result.size() - 1).getId())
        }
        return null
    }

    // Remapped variants due to Novel reference alleles (nra) can be deprecated
    void remediateDuplicateRemappedSSNotGeneratedByRemapping(List<? extends SubmittedVariantEntity> sves,
                                                             Class<? extends  SubmittedVariantEntity> collectionClassToUse) {
        List<Long> ssIDsToFind = sves.collect{sve -> sve.getAccession()}.unique()
        Query queryToGetSVEsAlongWithTheirDuplicates =
                query(this.criteriaToFindAssembly)
                        .addCriteria(this.criteriaToFindRemappedEntries)
                        .addCriteria(where("accession").in(ssIDsToFind))
        def svesAlongWithDuplicates =
                this.devMongoTemplate.find(queryToGetSVEsAlongWithTheirDuplicates, collectionClassToUse)
        // Source assembly here refers to the assembly that the SubmittedVariantEntity was remapped from
        svesAlongWithDuplicates.groupBy {sve -> sve.getRemappedFrom()}.each{
            sourceAssembly, svesWithDuplicates ->
                List<Long> ssIDsToFindInSourceAssembly = svesWithDuplicates.collect{sve -> sve.getAccession()}.unique()
                /** Remapped attributes collection looks like this
                 * {
                 "_id" : "GCA_000002035.2_996778337_CM002885.2_7587532_7587532_C_G",
                 "accession" : NumberLong(996778337),
                 "seq" : "GCA_000002035.2", (remapped)
                 "ref" : "C", (remapped)
                 "alt" : "G",(remapped)
                 "contig" : "CM002885.2", (remapped)
                 "start" : NumberLong(7587532), (remapped)
                 "end" : NumberLong(7587532), (remapped)
                 "type" : "SNV", (remapped)
                 "RS" : "rs513698182",
                 "st" : "-",
                 "CREATED" : "2014-05-23T10:30",
                 "PROJECT" : "NHGRI_DGS_NHGRI-1_FOUNDERS",
                 "AC" : "1",
                 "FILTER" : "PASS",
                 "rac" : "A-C",
                 "AF" : "1",
                 "TAX" : "7955",
                 "nra" : ""
                 }
                 */
                Query queryToFindSSIDsInSourceAssembly =
                        query(where("seq").is(sourceAssembly).and("accession").in(ssIDsToFindInSourceAssembly))
                def sourceAssemblyRemappingAttributesForSVEs =
                        this.devMongoTemplate.find(queryToFindSSIDsInSourceAssembly, Document.class,
                                this.remappedAttributesCollection)
                                .groupBy{doc -> "" + doc.get("seq") + "_" + doc.get("accession")}

                def sourceAssemblySVEEntries =
                        this.devMongoTemplate.find(queryToFindSSIDsInSourceAssembly.addCriteria(where("allelesMatch")
                                .exists(false).and("mapWeight").exists(false)), collectionClassToUse)
                                .groupBy{sve -> "" + sve.getReferenceSequenceAccession() + "_" + sve.getAccession()}

                def svesToDeprecate = svesWithDuplicates.groupBy{sve -> sve.getAccession()}
                        .collect{accession, svesWithTheAccession ->
                            def sourceAssemblyAndAccession = "" + "${sourceAssembly}_${accession}"
                            def correspondingRemappingAttributes =
                                    sourceAssemblyRemappingAttributesForSVEs.get(sourceAssemblyAndAccession)
                                            .collect{doc -> RemappingOutputParams.parseFrom(doc)}
                                            .unique()
                            def correspondingEntriesInSourceAssembly =
                                    sourceAssemblySVEEntries.get(sourceAssemblyAndAccession)
                            return determineSSWithInvalidRemaps(svesWithTheAccession, correspondingRemappingAttributes, correspondingEntriesInSourceAssembly)
                        }.flatten()
                this.ssDeprecationWriter.write(svesToDeprecate)
        }
    }

    // See https://docs.google.com/spreadsheets/d/1CW2V2uhsGDlTrDnu_NDpYgLUF-9LBSemDF61A7H7kYg/edit#gid=0
    static List<? extends SubmittedVariantEntity> determineSSWithInvalidRemaps(svesWithTheAccession, correspondingRemappingAttributes, correspondingEntriesInSourceAssembly) {
        if(svesWithTheAccession.size() == 0) return svesWithTheAccession
        def accession = ((SubmittedVariantEntity) (svesWithTheAccession.first())).getAccession()
        def assemblyToRemediate = ((SubmittedVariantEntity) (svesWithTheAccession.first())).getReferenceSequenceAccession()
        def svesWithTheAccessionGroupedByRefAlt = svesWithTheAccession.groupBy {sve ->
            def (normalizedRef, normalizedAlt) = normalizeAlleles(sve.getReferenceAllele(), sve.getAlternateAllele())
            return "${normalizedRef}_${normalizedAlt}"
        }
        def refAltPairsFromRemapping = correspondingEntriesInSourceAssembly.collect{sveFromSourceAssembly ->
            getPossibleRemappedAllelesForSVE(sveFromSourceAssembly, correspondingRemappingAttributes)}
                .collectMany{refAltPairList -> refAltPairList} // we need this hack to flatten "one level"
                .collect{refAltPair -> "${refAltPair[0]}_${refAltPair[1]}"}
        if (refAltPairsFromRemapping.size() == 0) {
            logger.error("No Ref/Alt remapping could be derived for accession ss$accession. This could be because this accession was not in the remapped VCF (and hence remapping attributes could not be determined) or remapping resulted in a non-variant....")
        }
        def svesNotGeneratedByRemapping =  svesWithTheAccessionGroupedByRefAlt.findAll{refAlt, sves ->
            // If the SVEs in the remapped assembly were not generated by remapping, deprecate them
            return !refAltPairsFromRemapping.contains(refAlt)
        }.values().flatten()
        if (svesNotGeneratedByRemapping.size() == svesWithTheAccessionGroupedByRefAlt.size()) {
            def remappedSVEWithSameTypeAsSourceSVE =
                    guessRemappedVariantWithSameTypeForSVE(svesWithTheAccession, correspondingEntriesInSourceAssembly)
            svesNotGeneratedByRemapping = svesWithTheAccession.findAll{sve -> !sve.equals(remappedSVEWithSameTypeAsSourceSVE)}
        }
        // If we can't trace the source variant for a remapped variant even after the above guess
        // report such a variant
        if (svesNotGeneratedByRemapping.size() == svesWithTheAccessionGroupedByRefAlt.size()) {
            logger.error("None of the submitted variants with accession ${accession} in assembly ${assemblyToRemediate} " +
                    "could be derived via remapping...")
        }
        return svesNotGeneratedByRemapping
    }

    // Transform a SVE based on the different remapping attributes and get the collection
    static List<String[]> getPossibleRemappedAllelesForSVE(sveFromSourceAssembly, correspondingRemappingAttributes) {
        return correspondingRemappingAttributes.collect {remappingAttribute ->
            def (refFromSourceAssembly, altFromSourceAssembly) =
            [sveFromSourceAssembly.getReferenceAllele(), sveFromSourceAssembly.getAlternateAllele()]
            // No transformations necessary when there is no novel reference or reference allele change or strand change
            if (Objects.isNull(remappingAttribute.referenceAlleleChange) && !remappingAttribute.novelReferenceAllele
                    && remappingAttribute.strand == "+") {
                return normalizeAlleles(refFromSourceAssembly, altFromSourceAssembly)
            }
            def (transformedRef, transformedAlt) = [refFromSourceAssembly, altFromSourceAssembly]
            if (remappingAttribute.strand == "-") {
                (transformedRef, transformedAlt) = [SequenceUtil.reverseComplement(transformedRef),
                                                    SequenceUtil.reverseComplement(transformedAlt)]
            }
            if (Objects.nonNull(remappingAttribute.referenceAlleleChange)) {
                def (referenceAlleleChangeOriginalRef, referenceAlleleChangeNewRef) = remappingAttribute.referenceAlleleChange.split("-")
                (transformedRef, transformedAlt) = applyReferenceAlleleChange(transformedRef, transformedAlt,
                        referenceAlleleChangeOriginalRef, referenceAlleleChangeNewRef)
            }
            if (transformedRef == transformedAlt) return null
            return normalizeAlleles(transformedRef, transformedAlt)
        }.findAll{refAltPair -> Objects.nonNull(refAltPair)}
    }

    static def applyReferenceAlleleChange(refFromSourceAssembly, altFromSourceAssembly,
                                          referenceAlleleChangeOriginalRef, referenceAlleleChangeNewRef) {
        // Try to check if there were "context" nucleotides added before or after the original reference from the source assembly
        def matcher = (referenceAlleleChangeOriginalRef =~ /(?<allelesToTheLeft>[A-Za-z]*)${refFromSourceAssembly}(?<allelesToTheRight>[A-Za-z]*)/)
        if (matcher.matches()) {
            def (allelesToTheLeft, allelesToTheRight) = [matcher.group("allelesToTheLeft"),
                                                         matcher.group("allelesToTheRight")]
            refFromSourceAssembly = "${allelesToTheLeft}${refFromSourceAssembly}${allelesToTheRight}"
            altFromSourceAssembly = "${allelesToTheLeft}${altFromSourceAssembly}${allelesToTheRight}"
        }

        refFromSourceAssembly = refFromSourceAssembly.replace(referenceAlleleChangeOriginalRef, referenceAlleleChangeNewRef)
        return [refFromSourceAssembly, altFromSourceAssembly]
    }

    static List<String> normalizeAlleles(String ref, String alt) {
        // use VariantCoreFields as an indirect way to normalize alleles
        def dummyVariant = new VariantCoreFields("dummyChr", 100L, ref, alt)
        return [dummyVariant.getReference(), dummyVariant.getAlternate()]
    }

    static SubmittedVariantEntity guessRemappedVariantWithSameTypeForSVE(List<? extends SubmittedVariantEntity> remappedSVEsWithTheAccession,
                                                                  List<? extends SubmittedVariantEntity> correspondingEntriesInSourceAssembly) {
        // This only works in a very specific case where
        // there is one variant in source assembly mapped to two or more variants in the remapped
        // and only one of the remapped variants has the same type as the variant in the source assembly
        // See Case 2 - https://docs.google.com/spreadsheets/d/1MTB1nwn7uooTk0uQzba_bmj5PTE6L_hhjKSj0Amj8Qo/edit#rangeid=1027635361
        def sourceVariantSVE = correspondingEntriesInSourceAssembly[0]
        if (correspondingEntriesInSourceAssembly.size() == 1) {
            VariantType sourceVariantType = convertSVEToVariant(sourceVariantSVE).getType()
            def remappedSVEsWithTheSameVariantType =
                    remappedSVEsWithTheAccession.findAll{sve -> isCompatibleRemappedVariantType(sourceVariantType,
                            convertSVEToVariant(sve).getType())}
            if (remappedSVEsWithTheSameVariantType.size() > 1) {
                logger.error("More than one remapped variant matched the variant type of source variant " +
                        "for accession ${sourceVariantSVE.getAccession()} " +
                        "in assembly ${sourceVariantSVE.getReferenceSequenceAccession()} when matching by type!!")
                return null
            }
            if (remappedSVEsWithTheSameVariantType.size() == 0) {
                logger.error("None of the remapped variants matched the variant type of the source variant " +
                        "for accession ${sourceVariantSVE.getAccession()} " +
                        "in assembly ${sourceVariantSVE.getReferenceSequenceAccession()} when matching by type!!")
                return null
            }
            return remappedSVEsWithTheSameVariantType[0]
        }
        else {
            logger.error("More than one source variant found for accession ${sourceVariantSVE.getAccession()} " +
                    "in assembly ${sourceVariantSVE.getReferenceSequenceAccession()} when attempting to match by type!!")
        }
        return null
    }

    static boolean isCompatibleRemappedVariantType(VariantType sourceVariantType, VariantType remappedVariantType) {
        if (sourceVariantType == remappedVariantType) {
            return true
        }
        if ((sourceVariantType.equals(VariantType.INS) && remappedVariantType.equals(VariantType.INDEL)) ||
                (sourceVariantType.equals(VariantType.DEL) && remappedVariantType.equals(VariantType.INDEL))) {
            return true
        }
        return false
    }

    static Variant convertSVEToVariant(SubmittedVariantEntity sve) {
        VariantCoreFields variantCoreFields = new VariantCoreFields(sve.getContig(), sve.getStart(),
                sve.getReferenceAllele(), sve.getAlternateAllele())
        return new Variant(variantCoreFields.getChromosome(), variantCoreFields.getStart(),
                variantCoreFields.getEnd(), variantCoreFields.getReference(), variantCoreFields.getAlternate())
    }
}


@EqualsAndHashCode
class RemappingOutputParams implements Comparable<RemappingOutputParams> {
    String seq
    Long accession
    String strand
    String referenceAlleleChange
    boolean novelReferenceAllele

    RemappingOutputParams() {

    }

    RemappingOutputParams(String seq, Long accession, String strand, String referenceAlleleChange,
                          boolean novelReferenceAllele) {
        this.seq = seq
        this.accession = accession
        this.strand = strand
        this.referenceAlleleChange = referenceAlleleChange
        this.novelReferenceAllele = novelReferenceAllele
    }

    static RemappingOutputParams parseFrom(Document remappingPropertiesDocument) {
        return new RemappingOutputParams(remappingPropertiesDocument.getString("seq"),
                remappingPropertiesDocument.getLong("accession"),
                remappingPropertiesDocument.get("st", "+").toString(),
                remappingPropertiesDocument.getString("rac"),
                remappingPropertiesDocument.containsKey("nra"))
    }

    @Override
    int compareTo(RemappingOutputParams other) {
        if(ObjectUtils.compare(this.seq, other.seq) == 0 &&
                ObjectUtils.compare(this.accession, other.accession) == 0 &&
                ObjectUtils.compare(this.strand, other.strand)  == 0 &&
                ObjectUtils.compare(this.referenceAlleleChange, other.referenceAlleleChange)  == 0 &&
                ObjectUtils.compare(this.novelReferenceAllele, other.novelReferenceAllele) == 0) {
            return 0
        }
        return -1
    }
}


@SpringBootApplication
class RemediateDuplicateRemappedSSApp  {

    static void main(String[] args) {
        ConfigurableApplicationContext context = SpringApplication.run(RemediateDuplicateRemappedSSApp.class, args)
        System.exit(SpringApplication.exit(context))
    }

}
