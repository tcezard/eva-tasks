package eva3399

import org.springframework.data.annotation.Id
import org.springframework.data.mongodb.core.mapping.Document
import uk.ac.ebi.eva.accession.core.model.eva.SubmittedVariantEntity

@Document
class ClashingSSHashes {
    @Id
    String ssHash
    List<SubmittedVariantEntity> clashingSS

    ClashingSSHashes() {}

    ClashingSSHashes(List<SubmittedVariantEntity> clashingSS) {
        this.clashingSS = clashingSS
        this.ssHash = this.clashingSS.get(0).hashedMessage
    }
}
