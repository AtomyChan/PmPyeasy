/*========================================================*/
/*                                                        */
/*  fitshedit.c     2005.12.19      version 4.1.1         */
/*                                                        */
/*  Copyright (C) 2005 by Wojtek Pych, CAMK PAN           */
/*                                                        */
/*  Writen for GNU project C and C++ Compiler.            */
/*                                                        */
/*  View and edit a FITS header with possible extensions. */
/*                                                        */
/*========================================================*/

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "errmess.h"
#include "pfitshead.h"

void usage(void);

/*------------------------------------------------------*/
void usage(void)
{
  printf("\n\t USAGE: fitshedit fits_file [keyword [value]]\n");
  exit(1);
}
/*------------------------------------------------------*/
char *test_keyword(char **keyword)
{
        char  *new_keyword;
        int   i;

  new_keyword=*keyword;

  if (strlen(new_keyword) > KEYWORD_SIZE)
  {
    printf("ERROR! %s : keyword too long\n", new_keyword);
    exit(2);
  }

  for (i=0; i<strlen(new_keyword); i++) if (islower(new_keyword[i]))
  {
    printf("ERROR! %s : lower case characters in the keyword\n", new_keyword);
    exit(3);
  }

  if (!(new_keyword=(char *)realloc(new_keyword, KEYWORD_SIZE*sizeof(char))))
    errmess("realloc(new_keyword)");

  for (i=strlen(new_keyword); i<KEYWORD_SIZE; i++) new_keyword[i]=' ';

  return(new_keyword);
}
/*------------------------------------------------------*/
char *test_val(char **value)
{
        char  *val,
              *new_val,
              isstr;        // 0: value=number  1: value=string
        int   i;

  val=*value;
  isstr=0;

  if (!isdigit(val[0]) && (val[0] != '+') && (val[0] != '-') && (val[0] != '.'))
    isstr=1;
  else
    for (i=1; i<strlen(val); i++)
      if (!isdigit(val[i]) && (val[i] != '.'))
        isstr=1;

  if (isstr)
  {
    if (!(new_val=(char *)calloc(strlen(val)+3, sizeof(char))))
      errmess("calloc(new_val)");
    new_val[0]='\'';
    memcpy((new_val+1), val, strlen(val));
    new_val[strlen(new_val)]='\'';
    new_val[strlen(new_val)]='\0';
  }
  else new_val=val;

  if (strlen(new_val) > VALUE_SIZE)
  {
    printf("ERROR! string:\n%s\nto long for FITS header card\n", new_val);
    exit(4);
  }

  return(new_val);
}
/*------------------------------------------------------*/
/*--------------------------------------------------------*/
/*  Return card number (starting from 0) or -1 on error   */
/*--------------------------------------------------------*/
int find_FITS_key(int hsize, char **header, int offset, char *keyword)
{
        int   i;              /*  loop numerator          */

  for (i=offset; i<hsize; i++)
    if (!strncmp(header[i], keyword, KEYWORD_SIZE)) break;

  if (i == hsize) return(-1);

  return(i);
}
/*------------------------------------------------------*/
char *get_comment(char *card)
{
        char  *comment;
        int   pc,
              j,
              len;

  pc=0;
  for (j=10; j<CARD_SIZE; j++)
  {
    if (card[j] == '/')
    {
      pc=j;
      break;
    }
  }

  if (!pc) return(NULL);

  len=CARD_SIZE-pc;
  if (!(comment=(char *)calloc(len, sizeof(char))))
    errmess("calloc(comment)");
  memcpy(comment, (card+pc), len);
  comment[len]='\0';

  return(comment);
}
/*------------------------------------------------------*/
int add_card(int *hsize, char ***header)
{
        char  **lheader,      // local pointer to the header
              card[CARD_SIZE];  // header card
        int   p,              // pointer to the card
              i,              // loop numerator
              lhsize;         // local header size

  lhsize=*hsize;
  lheader=*header;
  
  if ((p=get_FITS_key(lhsize, lheader, "END", card)) == -1)
  {
    printf("ERROR! END not found\n");
    exit(5);
  }

  memset(card, ' ', CARD_SIZE);
  if (!strncmp(card, lheader[p-1], CARD_SIZE))
  {
    for (i=0; i<p; i++)
      if (!strncmp(card, lheader[i], CARD_SIZE)) return(p);
  }

  if (!((p+1)%RECORD_CARDS))
  {
    {
      lhsize+=RECORD_CARDS;
      if (!(lheader=(char **)realloc(lheader, lhsize*sizeof(char *))))
        errmess("realloc(lheader)");
      for (i=1; i<=RECORD_CARDS; i++)
      {
        if (!(lheader[lhsize-i]=(char *)calloc(CARD_SIZE, sizeof(char))))
          errmess("malloc(lheader[lhsize-1])");
        memset(lheader[lhsize-i], ' ', CARD_SIZE);
      }
    }
  }

  memcpy(lheader[p+1], lheader[p], CARD_SIZE);
  memset(lheader[p], ' ', CARD_SIZE);

  *hsize=lhsize;
  *header=lheader;

  return(p);
}
/*------------------------------------------------------*/
char *edit_new_card(char *old_card, char *keyword, char *new_value)
{
        char  *new_card,      // new header card 
              *comment;       // card comment
        int   len,            // length of the card string
              clen;           // length of the comment

  if (!(new_card=(char *)calloc(CARD_SIZE, sizeof(char))))
    errmess("calloc(new_card)");
  memset(new_card, ' ', CARD_SIZE);

  sprintf(new_card, "%8s= ", keyword);

  len=strlen(new_value);
  memcpy((new_card+10), new_value, len);
  len+=10;
  if (len < CARD_SIZE)
  {
    new_card[len]=' ';
    len++;
  }

  if ((comment=get_comment(old_card)))
  {
    clen=strlen(comment);
    clen=((clen > CARD_SIZE-len) ? CARD_SIZE-len : clen);
    memcpy((new_card+len), comment, clen);
  }

  return(new_card);
}
/*------------------------------------------------------*/
int main(int argc, char *argv[])
{
        char  *name,          // input file name
              *keyword,       // keyword
              value[CARD_SIZE], // card valuer
              *new_card,      // new header card
              *new_value,     // new_value string
              ***header;      // header buffer
        int   i, j,           // loop numerators
              p,              // pointer to the card
              ne,             // number of extensions
              *hsize;         // vector of headers size
        long  *offsets;
        FILE  *inf;           // input file

  name = keyword = new_value = NULL;

/*** read command line ***/

  switch(argc)
  {
    case 4:
      if (!(new_value=(char *)calloc(strlen(argv[3])+1, sizeof(char))))
        errmess("calloc(new_value)");
      strcpy(new_value, argv[3]);
      new_value=test_val(&new_value);

    case 3:
      if (!(keyword=(char *)calloc(strlen(argv[2])+1, sizeof(char))))
        errmess("calloc(keyword)");
      strcpy(keyword, argv[2]);
      keyword=test_keyword(&keyword);

    case 2:
      name=argv[1];
      break;

    default: usage();
  }

/** read input file headers **/
  if (!(inf=fopen(name, "r"))) errmess(name);
  ne=read_FITS_headers(inf, &hsize, &header, &offsets);
  fclose(inf);

/** process headers **/
  switch(argc)
  {
    case 2:
      for (i=0; i<=ne; i++)
        for (j=0; j<hsize[i]; j++)
          printf("%s", header[i][j]);
      break;

    case 3:
      for (i=0; i<=ne; i++)
      {
        for (p=0; p!=-1; p++)
        {
          if ((p=find_FITS_key(hsize[i], header[i], p, keyword)) == -1) break;
          printf("%s", header[i][p]);
        }
        printf("\n");
      }
      break;

    case 4:
      if (!(inf=fopen(name, "r+"))) errmess(name);

      for (i=0; i<=ne; i++)
      {
        p=get_FITS_key(hsize[i], header[i], keyword, value);
        if ((p == -1)
            || !strcmp(keyword, "HISTORY ") || !strcmp(keyword, "COMMENT "))
          p=add_card(&hsize[i], &header[i]);

        if (!((p+1)%RECORD_CARDS))
        {
          printf("\tERROR! Must rewrite the whole file - not implemented.\n");
// Look into the trimhead.c to begin implementation.
          return(-2);
        }

        new_card=edit_new_card(header[i][p], keyword, new_value);
        memcpy(header[i][p], new_card, CARD_SIZE);

        printf("\n%s\n", header[i][p]);

        fseek(inf, offsets[i], SEEK_CUR);
        for (j=0; j<hsize[i]; j++)
          fwrite(header[i][j], sizeof(char), CARD_SIZE, inf);

        free(new_card);
      }

      printf("%s : header updated.\n", name);
      fclose(inf);
      break;

    default: printf("This should not happen.\n");
  }

  free(keyword);
  free(new_value);
  for (i=0; i<ne; i++)
  {
    for (j=0; j<hsize[i]; j++) free(header[i][j]);
    free(header[i]);
  }
  free(header);
  free(hsize);
  free(offsets);

  return(0);
}
/*** END ***/
