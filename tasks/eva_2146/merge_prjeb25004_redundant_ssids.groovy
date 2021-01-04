@Grab(group='uk.ac.ebi.eva', module='eva-accession-pipeline', version='0.4.5-SNAPSHOT')
@Grab('com.xlson.groovycsv:groovycsv:1.3')
import static com.xlson.groovycsv.CsvParser.parseCsv

import uk.ac.ebi.ampt2d.commons.accession.core.exceptions.AccessionMergedException
import uk.ac.ebi.eva.accession.core.configuration.nonhuman.MongoConfiguration
import uk.ac.ebi.eva.accession.core.configuration.nonhuman.SubmittedVariantAccessioningConfiguration
import uk.ac.ebi.eva.accession.core.configuration.ApplicationProperties
import uk.ac.ebi.eva.accession.core.service.nonhuman.SubmittedVariantAccessioningService


@Component
@Import(value=[SubmittedVariantAccessioningConfiguration.class, MongoConfiguration.class])
class MainApp implements CommandLineRunner {

	@Autowired
	private SubmittedVariantAccessioningService svAccessioningService;

	void run(String... args) {
		def file = args[0]
		def fileWithVariantsToMerge = new File(file).getText()

		def data = parseCsv(fileWithVariantsToMerge, separator: '\t', readFirstLine: true)
		data.each {line ->
			String fromSS = line[1]
			String toSS = line[5]
			merge(Long.parseLong(fromSS.replace("ss", "")), Long.parseLong(toSS.replace("ss", "")))
		}
	}

	void merge(Long fromSS, Long toSS) {
		try {
			svAccessioningService.merge(fromSS, toSS, "Identical submitted variant received multiple SS identifiers")
		}
		catch(AccessionMergedException ex) {
			println "$fromSS has already been merged to $toSS"
		}
	}
}
