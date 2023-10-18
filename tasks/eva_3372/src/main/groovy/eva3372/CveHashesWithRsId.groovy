package eva3372

import org.springframework.data.annotation.Id
import org.springframework.data.mongodb.core.mapping.Document

@Document
class CveHashesWithRsId {
    @Id
    String id

    String assembly

    Long rsId

    Set<String> cveHashes

    CveHashesWithRsId() {}

    CveHashesWithRsId(String id, String assembly, Long rsId, Set<String> cveHashes) {
        this.id = id
        this.assembly = assembly
        this.rsId = rsId
        this.cveHashes = cveHashes
    }
}
