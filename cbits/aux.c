#include "htslib/sam.h"

BGZF* get_bgzf(htsFile *h) { return h->fp.bgzf; }

char* get_header_text(bam_hdr_t *hdr) {
    return hdr->text;
}

uint32_t get_header_size(bam_hdr_t *hdr) {
    return hdr->l_text;
}

char* bam_chr(bam_hdr_t *h, int32_t i) {
    return h->target_name[i];
}

int bam_is_rev_(bam1_t *b) {
    return b->core.flag&BAM_FREVERSE != 0;
}

void bam_get_seq_(bam1_t *b, char *res, uint32_t l) {
    int32_t i;
    uint8_t *s = bam_get_seq(b);
    for (i = 0; i < l; ++i)
        res[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(s, i)];
}

int bam_get_qual_(bam1_t *b, char *res, uint32_t l) {
    int32_t i;
    uint8_t *s = bam_get_qual(b);
    if (s[0] == 0xff) return 1;
    for (i = 0; i < l; ++i) { res[i] = s[i]; }
    return 0;
}

void bam_get_cigar_(bam1_t *b, int *num, char *str, uint32_t l) {
    uint16_t i;
    uint32_t *cigar = bam_get_cigar(b);
    for (i = 0; i < l; ++i) {
        num[i] = bam_cigar_oplen(cigar[i]);
        str[i] = bam_cigar_opchr(cigar[i]);
    }
}

uint8_t* bam_get_aux_(bam1_t *b) { return bam_get_aux(b); }

int bam_get_l_aux_(bam1_t *b) { return bam_get_l_aux(b); }

