package eva3372

import org.springframework.data.annotation.Id
import org.springframework.data.mongodb.core.mapping.Document

@Document
class RsIdsWithCveHash {
    @Id
    String id

    String assembly

    String cveHash

    Set<Long> rsIds

    RsIdsWithCveHash() {}

    RsIdsWithCveHash(String id, String assembly, String cveHash, Set<Long> rsIds) {
        this.id = id
        this.assembly = assembly
        this.cveHash = cveHash
        this.rsIds = rsIds
    }
}
