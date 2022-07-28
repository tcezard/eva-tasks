package src

import org.springframework.data.annotation.Id
import org.springframework.data.mongodb.core.index.Indexed
import org.springframework.data.mongodb.core.mapping.Document
import uk.ac.ebi.ampt2d.commons.accession.hashing.SHA1HashingFunction

@Document
class DuplicateSSCategory {
    @Id
    String id

    @Indexed
    Long accession

    @Indexed
    String seq
    String category

    protected DuplicateSSCategory() {
    }

    DuplicateSSCategory(Long accession, String seq, String category) {
        this.accession = accession
        this.seq = seq
        this.category = category
        this.id = new SHA1HashingFunction().apply("${this.accession}_${this.seq}_${this.category}")
    }


    @Override
    public String toString() {
        return "DuplicateSSCategory{" +
                "id='" + id + '\'' +
                ", accession=" + accession +
                ", seq='" + seq + '\'' +
                ", category='" + category + '\'' +
                '}'
    }
}
