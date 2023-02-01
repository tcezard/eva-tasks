package eva3145

import org.springframework.data.annotation.Id
import org.springframework.data.mongodb.core.mapping.Document

@Document
class AsmMatchMetrics {
    @Id
    String assembly
    org.bson.Document metrics

    AsmMatchMetrics() {}

    AsmMatchMetrics(String assembly, Map<String, Long> metrics) {
        this.assembly = assembly
        this.metrics = new org.bson.Document(metrics)
    }
}