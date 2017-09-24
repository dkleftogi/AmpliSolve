#ifndef UTILITIES_H_INCLUDED
#define UTILITIES_H_INCLUDED

#include <assert.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include "samtools/sam.h"
#include "DataStructures.h"

static cstate_t g_cstate_null = { -1, 0, 0, 0 };

int isInteger(char *s)
{
    int i;
    for(i=0;i<strlen(s);i++)
      if(!(s[i] >= '0' && s[i] <= '9'))
        return 1;
    return 0;
}

int hasInteger(char *s)
{
    int i;
    for(i=0;i<strlen(s);i++)
      if((s[i] >= '0' && s[i] <= '9'))
        return 1;
    return 0;
}

void subSlash(char *s)
{
	int i;
    int pos=-1;
	#ifdef _WIN32
	for(i=0;i<strlen(s);i++)
      if(s[i] == '/')
        s[i]='\\';
	#else
    for(i=0;i<strlen(s);i++)
      if(s[i] == '\\')
        s[i]="/";
	#endif
}

int lastSlash(char *s)
{
    int i;
    int pos=-1;
	#ifdef _WIN32
	for(i=0;i<strlen(s);i++)
      if(s[i] == '\\')
        pos=i;
	#else
    for(i=0;i<strlen(s);i++)
      if(s[i] == '/')
        pos=i;
	#endif
	return pos;
}

int isKeyword(char *s)
{
    int i;
    if(strncmp(s,"STDOUT",6)==0 && strlen(s)==6)
        return 1;
    return 0;
}

int FileExists(const char * filename)
{
    FILE *istream;
    if ((istream=fopen(filename,"r"))==NULL)
    {
        return 1;
    }
    else
    {
        fclose ( istream );
        return 0;
    }
}

// Adaptation of SAM implementation
int resolve_cigar_mine(bam_pileup1_t *p, uint32_t pos,cstate_t *s,int *skip)
{
#define _cop(c) ((c)&BAM_CIGAR_MASK)
#define _cln(c) ((c)>>BAM_CIGAR_SHIFT)
	bam1_t *b = p->b;
	bam1_core_t *c = &b->core;
	uint32_t *cigar = bam1_cigar(b);
	int k, is_head = 0;
	int control = 0;
	while(control==0)
	{
	// determine the current CIGAR operation
    // fprintf(stderr, "%s\tpos=%d\tend=%d\t(%d,%d,%d)\n", bam1_qname(b), pos, s->end, s->k, s->x, s->y);
	if (s->k == -1) { // never processed
		is_head = 1;
		if (c->n_cigar == 1) { // just one operation, save a loop
		  if (_cop(cigar[0]) == BAM_CMATCH || _cop(cigar[0]) == BAM_CEQUAL || _cop(cigar[0]) == BAM_CDIFF) s->k = 0, s->x = c->pos, s->y = 0;
		} else { // find the first match or deletion
			for (k = 0, s->x = c->pos, s->y = 0; k < c->n_cigar; ++k) {
				int op = _cop(cigar[k]);
				int l = _cln(cigar[k]);
				if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CEQUAL || op == BAM_CDIFF) break;
				else if (op == BAM_CREF_SKIP) s->x += l;
				else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) s->y += l;
			}
			if(k>=c->n_cigar)
			{
                printf("%d - %d",k,c->n_cigar);
			}
			assert(k < c->n_cigar);
			s->k = k;
		}
	} else { // the read has been processed before
		int op, l = _cln(cigar[s->k]);
		if (pos - s->x >= l) { // jump to the next operation
			assert(s->k < c->n_cigar); // otherwise a bug: this function should not be called in this case
			op = _cop(cigar[s->k+1]);
			if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP || op == BAM_CEQUAL || op == BAM_CDIFF) { // jump to the next without a loop
			  if (_cop(cigar[s->k]) == BAM_CMATCH|| _cop(cigar[s->k]) == BAM_CEQUAL || _cop(cigar[s->k]) == BAM_CDIFF) s->y += l;
				s->x += l;
				++s->k;
			} else { // find the next M/D/N/=/X
			  if (_cop(cigar[s->k]) == BAM_CMATCH|| _cop(cigar[s->k]) == BAM_CEQUAL || _cop(cigar[s->k]) == BAM_CDIFF) s->y += l;
				s->x += l;
				for (k = s->k + 1; k < c->n_cigar; ++k) {
					op = _cop(cigar[k]), l = _cln(cigar[k]);
					if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP || op == BAM_CEQUAL || op == BAM_CDIFF) break;
					else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) s->y += l;
				}
				s->k = k;
			}
			assert(s->k < c->n_cigar); // otherwise a bug
		} // else, do nothing
		else
		{
		    if (cigar[s->k]&0xf != 0)
		    {
                *skip = 1;
		    }
            control=1;
		}
	}
	}
	{ // collect pileup information
		int op, l;
		op = _cop(cigar[s->k]); l = _cln(cigar[s->k]);
		p->is_del = p->indel = p->is_refskip = 0;
		if (s->x + l - 1 == pos && s->k + 1 < c->n_cigar) { // peek the next operation
			int op2 = _cop(cigar[s->k+1]);
			int l2 = _cln(cigar[s->k+1]);
			if (op2 == BAM_CDEL) p->indel = -(int)l2;
			else if (op2 == BAM_CINS) p->indel = l2;
			else if (op2 == BAM_CPAD && s->k + 2 < c->n_cigar) { // no working for adjacent padding
				int l3 = 0;
				for (k = s->k + 2; k < c->n_cigar; ++k) {
					op2 = _cop(cigar[k]); l2 = _cln(cigar[k]);
					if (op2 == BAM_CINS) l3 += l2;
					else if (op2 == BAM_CDEL || op2 == BAM_CMATCH || op2 == BAM_CREF_SKIP || op2 == BAM_CEQUAL || op2 == BAM_CDIFF) break;
				}
				if (l3 > 0) p->indel = l3;
			}
		}
		if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
			p->qpos = s->y + (pos - s->x);
		} else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
			p->is_del = 1; p->qpos = s->y; // FIXME: distinguish D and N!!!!!
			p->is_refskip = (op == BAM_CREF_SKIP);
		} // cannot be other operations; otherwise a bug
		p->is_head = (pos == c->pos); p->is_tail = (pos == s->end);
	}

	return 1;
}


#endif // UTILITIES_H_INCLUDED
