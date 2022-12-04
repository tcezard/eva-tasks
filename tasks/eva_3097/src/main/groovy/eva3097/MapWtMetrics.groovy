package eva3097

import org.springframework.data.annotation.Id
import org.springframework.data.mongodb.core.mapping.Document

@Document
class MapWtMetrics {
    @Id
    String assembly
    org.bson.Document metrics

    MapWtMetrics() {}

    MapWtMetrics(String assembly, Map<String, Integer> metrics) {
        this.assembly = assembly
        this.metrics = new org.bson.Document(metrics)
    }
}