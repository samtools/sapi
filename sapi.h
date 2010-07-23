#ifndef SAPI_H
#define SAPI_H

/* SAPI = Sequence Alignment Programming Interface */

/* Version: 0.0.3 */

#include <stdint.h>

/*********************
 * SAM bitwise flags *
 *********************/

#define SAM_FPAIRED        1
#define SAM_FPROPER_PAIR   2
#define SAM_FUNMAP         4
#define SAM_FMUNMAP        8
#define SAM_FREVERSE      16
#define SAM_FMREVERSE     32
#define SAM_FREAD1        64
#define SAM_FREAD2       128
#define SAM_FSECONDARY   256
#define SAM_FQCFAIL      512
#define SAM_FDUP        1024

/************************
 * SAM CIGAR operations *
 ************************/

#define SAM_CMATCH      0
#define SAM_CINS        1
#define SAM_CDEL        2
#define SAM_CREF_SKIP   3
#define SAM_CSOFT_CLIP  4
#define SAM_CHARD_CLIP  5
#define SAM_CPAD        6

/******************************
 * File open flags and others *
 ******************************/

#define SA_O_READ       1
#define SA_O_WRITE      2
#define SA_O_SAM        0x10000
#define SA_O_BAM        0x20000
#define SA_O_BIOHDF     0x40000

#define SA_FREE_HEADER  0x1
#define SA_FREE_INDEX   0x2
#define SA_FREE_ALL     0xffffffffu

#define SA_GET_CORE     0x1
#define SA_GET_CIGAR    0x2
#define SA_GET_QNAME    0x4
#define SA_GET_SEQ      0x8
#define SA_GET_QUAL     0x10
#define SA_GET_TAGS     0x20
#define SA_GET_ALL      0xffffffffu

/*************
 * Alignment *
 *************/

typedef struct {
	struct {
		int32_t tid;
		int32_t pos;
		uint32_t dummay:16, qual:8, l_qname:8; // dummy for BAM
		uint32_t flag:16, n_cigar:16;
		int32_t l_qseq;
		int32_t mtid;
		int32_t mpos;
		int32_t isize;
	} core;
	int l_data, m_data;
	uint8_t *data;
} sa_aln_t;

#define sa_aln_strand(b) (((b)->core.flag&BAM_FREVERSE) != 0)
#define sa_aln_mstrand(b) (((b)->core.flag&BAM_FMREVERSE) != 0)
#define sa_aln_cigar(b) ((uint32_t*)((b)->data + (b)->core.l_qname))
#define sa_aln_qname(b) ((char*)((b)->data))
#define sa_aln_seq(b) ((b)->data + (b)->core.n_cigar*4 + (b)->core.l_qname)
#define sa_aln_seqi(s, i) ((s)[(i)/2] >> 4*(1-(i)%2) & 0xf)
#define sa_aln_qual(b) ((b)->data + (b)->core.n_cigar*4 + (b)->core.l_qname + ((b)->core.l_qseq + 1)/2)
#define sa_aln_aux(b) ((b)->data + (b)->core.n_cigar*4 + (b)->core.l_qname + (b)->core.l_qseq + ((b)->core.l_qseq + 1)/2)

/*************************
 * Format dependent APIs *
 *************************/

typedef void sa_file_t;
typedef void sa_itr_t;
typedef int (*sa_hook_f)(sa_aln_t *aln, void *data);

typedef struct {
	int n_ref, l_text;
	int *ref_len;
	char **ref_name;
	char *text;
} sa_hdrinfo_t;

#ifdef __cplusplus
extern "C" {
#endif

	/* Open an alignment file */
	sa_file_t *sa_open(const char *fn, int mode);

	/* Close the file */
	int sa_close(sa_file_t *fp);

	/* Free part of internal data in fp (to save memory) */
	int sa_free(sa_file_t *fp, int which);

	/* Get an iterator from the current file position */
	sa_itr_t *sa_query_current(sa_file_t *fp);

	/* Get an iterator for alignments starting from a coordinate */
	sa_itr_t *sa_query_start(sa_file_t *fp, int ref, int beg);

	/* Get an iterator for alignments overlapping a region */
	sa_itr_t *sa_query_overlap(sa_file_t *fp, int ref, int beg, int end);

	/* Get an iterator for alignments in an array */
	sa_itr_t *sa_query_array(int n, sa_aln_t *array);

	/* Set hook which is called when an alignment is read by the iterator */
	int sa_set_hook(sa_itr_t *itr, sa_hook_f func, void *data);

	/* What information to retrieve (for column-sorted format, retrieving all information is inefficient) */
	int sa_set_content(sa_itr_t *itr, int which);

	/* Destroy an iterator */
	int sa_itr_destroy(sa_itr_t *itr);

	/* Read the next alignment from the iterator */
	int sa_next(sa_itr_t *itr, sa_aln_t *aln);

	/* Write an alignment to a file */
	int sa_write(sa_file_t *fp, const sa_aln_t *aln);

	/* Write the header */
	int sa_write_header(sa_file_t *fp);

	/* Get key header information */
	const sa_hdrinfo_t *sa_get_hdrinfo(const sa_file_t *fp);

	/* Set header information */
	int sa_set_hdrinfo(sa_file_t *fp, const sa_hdrinfo_t *header);

	/* Convert string reference name to integer ID */
	int sa_refname2id(sa_file_t *fp, const char *refname);

	int sa_get_supported_index_types(sa_file_t *fp, /*OUT*/ int **indexes, /*OUT*/ int *n_indexes);
	int sa_build_index(sa_file_t *fp, int index_type);
	int sa_get_indexes(sa_file_t *fp, /*OUT*/ int **indexes, /*OUT*/ int *n_indexes);
	int sa_use_index(sa_file_t *fp, int index_type);
	int sa_unload_index(sa_file_t *fp);

#ifdef __cplusplus
}
#endif

/***************************
 * Format independent APIs *
 ***************************/

typedef void sa_plp_t;

typedef struct {
	sa_aln_t *b;
	int32_t qpos;
	int indel, level;
	uint32_t is_del:1, is_head:1, is_tail:1;
} sa_plpinfo_t;

#define SA_TTYPE_INT   1
#define SA_TTYPE_FLOAT 2
#define SA_TTYPE_CHAR  3
#define SA_TTYPE_STR   4

typedef struct {
	int type; // or use enum
	union {
		int i;
		char c;
		float f;
		const char *Z;
	} val;
} sa_tagval_t;

#ifdef __cplusplus
extern "C" {
#endif

	sa_aln_t *sa_aln_dup(const sa_aln_t *aln);
	void sa_aln_destroy(sa_aln_t *aln);

	/* Get the alignment in the SAM format */
	char *sa_aln2sam(const sa_aln_t *aln);

	/* Parse a SAM alignment line */
	int sa_sam2aln(const char *sam, sa_aln_t *aln);

	/* Initiate a pileup iterator */
	sa_plp_t *sa_plp_init(sa_itr_t *itr);

	/* Get the next pileup line; must not be modified by the caller */
	const sa_plpinfo_t *sa_plp_next(sa_plp_t *plp, int *n);

	/* Destroy a pileup iterator */
	void sa_plp_destroy(sa_plp_t *plp);

	int sa_get_tag(sa_aln_t *aln, const char tag[2], sa_tagval_t *tv);
	int sa_set_tag(sa_aln_t *aln, const char tag[2], const sa_tagval_t *tv);

#ifdef __cplusplus
}
#endif

#endif
