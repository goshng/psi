/* seq.c -- biological sequence
   Copyright (C) 2006 Sang Chul Choi
  
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.
  
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
  
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
   MA 02110-1301, USA.
*/

/** @start 1 */
#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "common.h"
#include "error.h"
#include "seq.h"

static int translation_table_id = PSI_CODON_TABLE_1;
static int init_codon = 0;

/* char-to-int value of amino acid like INT_AMINOACID[(int)'A'] => PSI_AA_ALA */
const int INT_AMINOACID[128] = {
    21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 
    21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 
  /*     !   "   #   $   %   &   '   (   )   *   +   ,   -   .   / */
    21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 20, 21, 21, 21, 21, 21, 
  /* 0   1   2   3   4   5   6   7   8   9   :   ;   <   =   >   ? */
    21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 
  /* @   A   B   C   D   E   F   G   H   I   J   K   L   M   N   O */
    21,  0, 21,  4,  3,  6, 13,  7,  8,  9, 21, 11, 10, 12,  2, 21, 
  /* P   Q   R   S   T   U   V   W   X   Y   Z   [   \   ]   ^   _ */
    14,  5,  1, 15, 16, 21, 19, 17, 21, 18, 21, 21, 21, 21, 21, 21, 
  /* `   a   b   c   d   e   f   g   h   i   j   k   l   m   n   o */
    21,  0, 21,  4,  3,  6, 13,  7,  8,  9, 21, 11, 10, 12,  2, 21, 
  /* p   q   r   s   t   u   v   w   x   y   z   {   |   }   ~   del */
    14,  5,  1, 15, 16, 21, 19, 17, 21, 18, 21, 21, 21, 21, 21, 21
};

/* char-to-int value of nucleotide like INT_NUCLEOTIDE[(int)'A'] => PSI_DNA_A */
const int INT_NUCLEOTIDE[128] = {
     4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  
     4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  
     4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  
     4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  
  /* @   A   B   C   D   E   F   G   H   I   J   K   L   M   N   O */
     4,  0,  4,  1,  4,  4,  4,  2,  4,  4,  4,  4,  4,  4,  4,  4,  
  /* P   Q   R   S   T   U   V   W   X   Y   Z   [   \   ]   ^   _ */
     4,  4,  4,  4,  3,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  
  /* `   a   b   c   d   e   f   g   h   i   j   k   l   m   n   o */
     4,  0,  4,  1,  4,  4,  4,  2,  4,  4,  4,  4,  4,  4,  4,  4,  
  /* p   q   r   s   t   u   v   w   x   y   z   {   |   }   ~   del */
     4,  4,  4,  4,  3,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4
};

/* int-to-char value of nucleotide like NUCLEOTIDE[PSI_DNA_A] ==> 'A' */
const char NUCLEOTIDE[4] = { 'A', 'C', 'G', 'T' };

/* int-to-char value of amino acid like AMINOACID[PSI_AA_ALA] ==> 'A' */
const char AMINOACID[21] = { 
   'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
/*  65,  82,  78,  68,  67,  81,  69,  71,  72,  73, */
/*  97, 114, 110, 100,  99, 113, 101, 103, 104, 105, */
   'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V',
/*  76,  75,  77,  70,  80,  83,  84,  87,  89,  86  */
/* 108, 107, 109, 102, 112, 115, 116, 119, 121, 118  */
   '*'
/*  42 */
};
 
/* The translation table 
   Nuc2AATable[PSI_DNA_A][PSI_DNA_A][PSI_DNA_A] ==> PSI_AA_LYS
   The following table is for the standard codon table 
   and we will read the "gc.prt" table to change this table.
   There were 17 codon tables available as of this writing.
   We will have a function to parse and populate this 3-D array
   with it. We will hard-code the information in a 4-D array.
*/
int Nuc2AATable[4][4][4] = 
{
  {
      { 11,  2, 11,  2 },
      { 16, 16, 16, 16 },
      {  1, 15,  1, 15 },
      {  9,  9, 12,  9 }
  },
  {
      {  5,  8,  5,  8 },
      { 14, 14, 14, 14 },
      {  1,  1,  1,  1 },
      { 10, 10, 10, 10 }
  },
  {
      {  6,  3,  6,  3 },
      {  0,  0,  0,  0 },
      {  7,  7,  7,  7 },
      { 19, 19, 19, 19 }
  },
  {
      { 20, 18, 20, 18 },
      { 15, 15, 15, 15 },
      { 20,  4, 17,  4 },
      { 10, 13, 10, 13 }
  }
};

int Nuc2AATableStart[4][4][4] = 
{
  {
    { 11,   2,  11,   2 },
    { 16,  16,  16,  16 },
    {  1,  15,   1,  15 },
    {  9,   9,  12,   9 }
  },
  {
    {  5,   8,   5,   8 },
    { 14,  14,  14,  14 },
    {  1,   1,   1,   1 },
    { 10,  10,  12,  10 }
  },
  {
    {  6,   3,   6,   3 },
    {  0,   0,   0,   0 },
    {  7,   7,   7,   7 },
    { 19,  19,  19,  19 }
  },
  {
    { 20,  18,  20,  18 },
    { 15,  15,  15,  15 },
    { 20,   4,  17,   4 },
    { 10,  13,  12,  13 }
  }
};

const int Nuc2AATables[17][4][4][4] = 
{
  {
    {
      { 11,   2,  11,   2 },
      { 16,  16,  16,  16 },
      {  1,  15,   1,  15 },
      {  9,   9,  12,   9 }
    },
    {
      {  5,   8,   5,   8 },
      { 14,  14,  14,  14 },
      {  1,   1,   1,   1 },
      { 10,  10,  10,  10 }
    },
    {
      {  6,   3,   6,   3 },
      {  0,   0,   0,   0 },
      {  7,   7,   7,   7 },
      { 19,  19,  19,  19 }
    },
    {
      { 20,  18,  20,  18 },
      { 15,  15,  15,  15 },
      { 20,   4,  17,   4 },
      { 10,  13,  10,  13 }
    }
  },
  {
    {
      { 11,   2,  11,   2 },
      { 16,  16,  16,  16 },
      { 20,  15,  20,  15 },
      { 12,   9,  12,   9 }
    },
    {
      {  5,   8,   5,   8 },
      { 14,  14,  14,  14 },
      {  1,   1,   1,   1 },
      { 10,  10,  10,  10 }
    },
    {
      {  6,   3,   6,   3 },
      {  0,   0,   0,   0 },
      {  7,   7,   7,   7 },
      { 19,  19,  19,  19 }
    },
    {
      { 20,  18,  20,  18 },
      { 15,  15,  15,  15 },
      { 17,   4,  17,   4 },
      { 10,  13,  10,  13 }
    }
  },
  {
    {
      { 11,   2,  11,   2 },
      { 16,  16,  16,  16 },
      {  1,  15,   1,  15 },
      { 12,   9,  12,   9 }
    },
    {
      {  5,   8,   5,   8 },
      { 14,  14,  14,  14 },
      {  1,   1,   1,   1 },
      { 16,  16,  16,  16 }
    },
    {
      {  6,   3,   6,   3 },
      {  0,   0,   0,   0 },
      {  7,   7,   7,   7 },
      { 19,  19,  19,  19 }
    },
    {
      { 20,  18,  20,  18 },
      { 15,  15,  15,  15 },
      { 17,   4,  17,   4 },
      { 10,  13,  10,  13 }
    }
  },
  {
    {
      { 11,   2,  11,   2 },
      { 16,  16,  16,  16 },
      {  1,  15,   1,  15 },
      {  9,   9,  12,   9 }
    },
    {
      {  5,   8,   5,   8 },
      { 14,  14,  14,  14 },
      {  1,   1,   1,   1 },
      { 10,  10,  10,  10 }
    },
    {
      {  6,   3,   6,   3 },
      {  0,   0,   0,   0 },
      {  7,   7,   7,   7 },
      { 19,  19,  19,  19 }
    },
    {
      { 20,  18,  20,  18 },
      { 15,  15,  15,  15 },
      { 17,   4,  17,   4 },
      { 10,  13,  10,  13 }
    }
  },
  {
    {
      { 11,   2,  11,   2 },
      { 16,  16,  16,  16 },
      { 15,  15,  15,  15 },
      { 12,   9,  12,   9 }
    },
    {
      {  5,   8,   5,   8 },
      { 14,  14,  14,  14 },
      {  1,   1,   1,   1 },
      { 10,  10,  10,  10 }
    },
    {
      {  6,   3,   6,   3 },
      {  0,   0,   0,   0 },
      {  7,   7,   7,   7 },
      { 19,  19,  19,  19 }
    },
    {
      { 20,  18,  20,  18 },
      { 15,  15,  15,  15 },
      { 17,   4,  17,   4 },
      { 10,  13,  10,  13 }
    }
  },
  {
    {
      { 11,   2,  11,   2 },
      { 16,  16,  16,  16 },
      {  1,  15,   1,  15 },
      {  9,   9,  12,   9 }
    },
    {
      {  5,   8,   5,   8 },
      { 14,  14,  14,  14 },
      {  1,   1,   1,   1 },
      { 10,  10,  10,  10 }
    },
    {
      {  6,   3,   6,   3 },
      {  0,   0,   0,   0 },
      {  7,   7,   7,   7 },
      { 19,  19,  19,  19 }
    },
    {
      {  5,  18,   5,  18 },
      { 15,  15,  15,  15 },
      { 20,   4,  17,   4 },
      { 10,  13,  10,  13 }
    }
  },
  {
    {
      {  2,   2,  11,   2 },
      { 16,  16,  16,  16 },
      { 15,  15,  15,  15 },
      {  9,   9,  12,   9 }
    },
    {
      {  5,   8,   5,   8 },
      { 14,  14,  14,  14 },
      {  1,   1,   1,   1 },
      { 10,  10,  10,  10 }
    },
    {
      {  6,   3,   6,   3 },
      {  0,   0,   0,   0 },
      {  7,   7,   7,   7 },
      { 19,  19,  19,  19 }
    },
    {
      { 20,  18,  20,  18 },
      { 15,  15,  15,  15 },
      { 17,   4,  17,   4 },
      { 10,  13,  10,  13 }
    }
  },
  {
    {
      { 11,   2,  11,   2 },
      { 16,  16,  16,  16 },
      {  1,  15,   1,  15 },
      {  9,   9,  12,   9 }
    },
    {
      {  5,   8,   5,   8 },
      { 14,  14,  14,  14 },
      {  1,   1,   1,   1 },
      { 10,  10,  10,  10 }
    },
    {
      {  6,   3,   6,   3 },
      {  0,   0,   0,   0 },
      {  7,   7,   7,   7 },
      { 19,  19,  19,  19 }
    },
    {
      { 20,  18,  20,  18 },
      { 15,  15,  15,  15 },
      {  4,   4,  17,   4 },
      { 10,  13,  10,  13 }
    }
  },
  {
    {
      { 11,   2,  11,   2 },
      { 16,  16,  16,  16 },
      {  1,  15,   1,  15 },
      {  9,   9,  12,   9 }
    },
    {
      {  5,   8,   5,   8 },
      { 14,  14,  14,  14 },
      {  1,   1,   1,   1 },
      { 10,  10,  10,  10 }
    },
    {
      {  6,   3,   6,   3 },
      {  0,   0,   0,   0 },
      {  7,   7,   7,   7 },
      { 19,  19,  19,  19 }
    },
    {
      { 20,  18,  20,  18 },
      { 15,  15,  15,  15 },
      { 20,   4,  17,   4 },
      { 10,  13,  10,  13 }
    }
  },
  {
    {
      { 11,   2,  11,   2 },
      { 16,  16,  16,  16 },
      {  1,  15,   1,  15 },
      {  9,   9,  12,   9 }
    },
    {
      {  5,   8,   5,   8 },
      { 14,  14,  14,  14 },
      {  1,   1,   1,   1 },
      { 10,  10,  15,  10 }
    },
    {
      {  6,   3,   6,   3 },
      {  0,   0,   0,   0 },
      {  7,   7,   7,   7 },
      { 19,  19,  19,  19 }
    },
    {
      { 20,  18,  20,  18 },
      { 15,  15,  15,  15 },
      { 20,   4,  17,   4 },
      { 10,  13,  10,  13 }
    }
  },
  {
    {
      { 11,   2,  11,   2 },
      { 16,  16,  16,  16 },
      {  7,  15,   7,  15 },
      { 12,   9,  12,   9 }
    },
    {
      {  5,   8,   5,   8 },
      { 14,  14,  14,  14 },
      {  1,   1,   1,   1 },
      { 10,  10,  10,  10 }
    },
    {
      {  6,   3,   6,   3 },
      {  0,   0,   0,   0 },
      {  7,   7,   7,   7 },
      { 19,  19,  19,  19 }
    },
    {
      { 20,  18,  20,  18 },
      { 15,  15,  15,  15 },
      { 17,   4,  17,   4 },
      { 10,  13,  10,  13 }
    }
  },
  {
    {
      {  2,   2,  11,   2 },
      { 16,  16,  16,  16 },
      { 15,  15,  15,  15 },
      {  9,   9,  12,   9 }
    },
    {
      {  5,   8,   5,   8 },
      { 14,  14,  14,  14 },
      {  1,   1,   1,   1 },
      { 10,  10,  10,  10 }
    },
    {
      {  6,   3,   6,   3 },
      {  0,   0,   0,   0 },
      {  7,   7,   7,   7 },
      { 19,  19,  19,  19 }
    },
    {
      { 18,  18,  20,  18 },
      { 15,  15,  15,  15 },
      { 17,   4,  17,   4 },
      { 10,  13,  10,  13 }
    }
  },
  {
    {
      { 11,   2,  11,   2 },
      { 16,  16,  16,  16 },
      {  1,  15,   1,  15 },
      {  9,   9,  12,   9 }
    },
    {
      {  5,   8,   5,   8 },
      { 14,  14,  14,  14 },
      {  1,   1,   1,   1 },
      { 10,  10,  10,  10 }
    },
    {
      {  6,   3,   6,   3 },
      {  0,   0,   0,   0 },
      {  7,   7,   7,   7 },
      { 19,  19,  19,  19 }
    },
    {
      { 20,  18,   5,  18 },
      { 15,  15,  15,  15 },
      { 20,   4,  17,   4 },
      { 10,  13,  10,  13 }
    }
  },
  {
    {
      { 11,   2,  11,   2 },
      { 16,  16,  16,  16 },
      {  1,  15,   1,  15 },
      {  9,   9,  12,   9 }
    },
    {
      {  5,   8,   5,   8 },
      { 14,  14,  14,  14 },
      {  1,   1,   1,   1 },
      { 10,  10,  10,  10 }
    },
    {
      {  6,   3,   6,   3 },
      {  0,   0,   0,   0 },
      {  7,   7,   7,   7 },
      { 19,  19,  19,  19 }
    },
    {
      { 20,  18,  10,  18 },
      { 15,  15,  15,  15 },
      { 20,   4,  17,   4 },
      { 10,  13,  10,  13 }
    }
  },
  {
    {
      {  2,   2,  11,   2 },
      { 16,  16,  16,  16 },
      { 15,  15,  15,  15 },
      { 12,   9,  12,   9 }
    },
    {
      {  5,   8,   5,   8 },
      { 14,  14,  14,  14 },
      {  1,   1,   1,   1 },
      { 10,  10,  10,  10 }
    },
    {
      {  6,   3,   6,   3 },
      {  0,   0,   0,   0 },
      {  7,   7,   7,   7 },
      { 19,  19,  19,  19 }
    },
    {
      { 20,  18,  20,  18 },
      { 15,  15,  15,  15 },
      { 17,   4,  17,   4 },
      { 10,  13,  10,  13 }
    }
  },
  {
    {
      { 11,   2,  11,   2 },
      { 16,  16,  16,  16 },
      {  1,  15,   1,  15 },
      {  9,   9,  12,   9 }
    },
    {
      {  5,   8,   5,   8 },
      { 14,  14,  14,  14 },
      {  1,   1,   1,   1 },
      { 10,  10,  10,  10 }
    },
    {
      {  6,   3,   6,   3 },
      {  0,   0,   0,   0 },
      {  7,   7,   7,   7 },
      { 19,  19,  19,  19 }
    },
    {
      { 20,  18,  10,  18 },
      { 20,  15,  15,  15 },
      { 20,   4,  17,   4 },
      { 10,  13,  10,  13 }
    }
  },
  {
    {
      { 11,   2,  11,   2 },
      { 16,  16,  16,  16 },
      {  1,  15,   1,  15 },
      {  9,   9,  12,   9 }
    },
    {
      {  5,   8,   5,   8 },
      { 14,  14,  14,  14 },
      {  1,   1,   1,   1 },
      { 10,  10,  10,  10 }
    },
    {
      {  6,   3,   6,   3 },
      {  0,   0,   0,   0 },
      {  7,   7,   7,   7 },
      { 19,  19,  19,  19 }
    },
    {
      { 20,  18,  20,  18 },
      { 15,  15,  15,  15 },
      { 20,   4,  17,   4 },
      { 20,  13,  10,  13 }
    }
  }
};

const int Nuc2AATablesStart[17][4][4][4] = 
{
  {
    {
      { 11,   2,  11,   2 },
      { 16,  16,  16,  16 },
      {  1,  15,   1,  15 },
      {  9,   9,  12,   9 }
    },
    {
      {  5,   8,   5,   8 },
      { 14,  14,  14,  14 },
      {  1,   1,   1,   1 },
      { 10,  10,  12,  10 }
    },
    {
      {  6,   3,   6,   3 },
      {  0,   0,   0,   0 },
      {  7,   7,   7,   7 },
      { 19,  19,  19,  19 }
    },
    {
      { 20,  18,  20,  18 },
      { 15,  15,  15,  15 },
      { 20,   4,  17,   4 },
      { 10,  13,  12,  13 }
    }
  },
  {
    {
      { 11,   2,  11,   2 },
      { 16,  16,  16,  16 },
      { 20,  15,  20,  15 },
      { 12,  12,  12,  12 }
    },
    {
      {  5,   8,   5,   8 },
      { 14,  14,  14,  14 },
      {  1,   1,   1,   1 },
      { 10,  10,  10,  10 }
    },
    {
      {  6,   3,   6,   3 },
      {  0,   0,   0,   0 },
      {  7,   7,   7,   7 },
      { 19,  19,  12,  19 }
    },
    {
      { 20,  18,  20,  18 },
      { 15,  15,  15,  15 },
      { 17,   4,  17,   4 },
      { 10,  13,  10,  13 }
    }
  },
  {
    {
      { 11,   2,  11,   2 },
      { 16,  16,  16,  16 },
      {  1,  15,   1,  15 },
      { 12,   9,  12,   9 }
    },
    {
      {  5,   8,   5,   8 },
      { 14,  14,  14,  14 },
      {  1,   1,   1,   1 },
      { 16,  16,  16,  16 }
    },
    {
      {  6,   3,   6,   3 },
      {  0,   0,   0,   0 },
      {  7,   7,   7,   7 },
      { 19,  19,  19,  19 }
    },
    {
      { 20,  18,  20,  18 },
      { 15,  15,  15,  15 },
      { 17,   4,  17,   4 },
      { 10,  13,  10,  13 }
    }
  },
  {
    {
      { 11,   2,  11,   2 },
      { 16,  16,  16,  16 },
      {  1,  15,   1,  15 },
      { 12,  12,  12,  12 }
    },
    {
      {  5,   8,   5,   8 },
      { 14,  14,  14,  14 },
      {  1,   1,   1,   1 },
      { 10,  10,  12,  10 }
    },
    {
      {  6,   3,   6,   3 },
      {  0,   0,   0,   0 },
      {  7,   7,   7,   7 },
      { 19,  19,  12,  19 }
    },
    {
      { 20,  18,  20,  18 },
      { 15,  15,  15,  15 },
      { 17,   4,  17,   4 },
      { 12,  13,  12,  13 }
    }
  },
  {
    {
      { 11,   2,  11,   2 },
      { 16,  16,  16,  16 },
      { 15,  15,  15,  15 },
      { 12,  12,  12,  12 }
    },
    {
      {  5,   8,   5,   8 },
      { 14,  14,  14,  14 },
      {  1,   1,   1,   1 },
      { 10,  10,  10,  10 }
    },
    {
      {  6,   3,   6,   3 },
      {  0,   0,   0,   0 },
      {  7,   7,   7,   7 },
      { 19,  19,  12,  19 }
    },
    {
      { 20,  18,  20,  18 },
      { 15,  15,  15,  15 },
      { 17,   4,  17,   4 },
      { 10,  13,  12,  13 }
    }
  },
  {
    {
      { 11,   2,  11,   2 },
      { 16,  16,  16,  16 },
      {  1,  15,   1,  15 },
      {  9,   9,  12,   9 }
    },
    {
      {  5,   8,   5,   8 },
      { 14,  14,  14,  14 },
      {  1,   1,   1,   1 },
      { 10,  10,  10,  10 }
    },
    {
      {  6,   3,   6,   3 },
      {  0,   0,   0,   0 },
      {  7,   7,   7,   7 },
      { 19,  19,  19,  19 }
    },
    {
      {  5,  18,   5,  18 },
      { 15,  15,  15,  15 },
      { 20,   4,  17,   4 },
      { 10,  13,  10,  13 }
    }
  },
  {
    {
      {  2,   2,  11,   2 },
      { 16,  16,  16,  16 },
      { 15,  15,  15,  15 },
      {  9,   9,  12,   9 }
    },
    {
      {  5,   8,   5,   8 },
      { 14,  14,  14,  14 },
      {  1,   1,   1,   1 },
      { 10,  10,  10,  10 }
    },
    {
      {  6,   3,   6,   3 },
      {  0,   0,   0,   0 },
      {  7,   7,   7,   7 },
      { 19,  19,  12,  19 }
    },
    {
      { 20,  18,  20,  18 },
      { 15,  15,  15,  15 },
      { 17,   4,  17,   4 },
      { 10,  13,  10,  13 }
    }
  },
  {
    {
      { 11,   2,  11,   2 },
      { 16,  16,  16,  16 },
      {  1,  15,   1,  15 },
      {  9,   9,  12,   9 }
    },
    {
      {  5,   8,   5,   8 },
      { 14,  14,  14,  14 },
      {  1,   1,   1,   1 },
      { 10,  10,  10,  10 }
    },
    {
      {  6,   3,   6,   3 },
      {  0,   0,   0,   0 },
      {  7,   7,   7,   7 },
      { 19,  19,  19,  19 }
    },
    {
      { 20,  18,  20,  18 },
      { 15,  15,  15,  15 },
      {  4,   4,  17,   4 },
      { 10,  13,  10,  13 }
    }
  },
  {
    {
      { 11,   2,  11,   2 },
      { 16,  16,  16,  16 },
      {  1,  15,   1,  15 },
      { 12,  12,  12,  12 }
    },
    {
      {  5,   8,   5,   8 },
      { 14,  14,  14,  14 },
      {  1,   1,   1,   1 },
      { 10,  10,  12,  10 }
    },
    {
      {  6,   3,   6,   3 },
      {  0,   0,   0,   0 },
      {  7,   7,   7,   7 },
      { 19,  19,  12,  19 }
    },
    {
      { 20,  18,  20,  18 },
      { 15,  15,  15,  15 },
      { 20,   4,  17,   4 },
      { 10,  13,  12,  13 }
    }
  },
  {
    {
      { 11,   2,  11,   2 },
      { 16,  16,  16,  16 },
      {  1,  15,   1,  15 },
      {  9,   9,  12,   9 }
    },
    {
      {  5,   8,   5,   8 },
      { 14,  14,  14,  14 },
      {  1,   1,   1,   1 },
      { 10,  10,  12,  10 }
    },
    {
      {  6,   3,   6,   3 },
      {  0,   0,   0,   0 },
      {  7,   7,   7,   7 },
      { 19,  19,  19,  19 }
    },
    {
      { 20,  18,  20,  18 },
      { 15,  15,  15,  15 },
      { 20,   4,  17,   4 },
      { 10,  13,  10,  13 }
    }
  },
  {
    {
      { 11,   2,  11,   2 },
      { 16,  16,  16,  16 },
      {  7,  15,   7,  15 },
      { 12,   9,  12,   9 }
    },
    {
      {  5,   8,   5,   8 },
      { 14,  14,  14,  14 },
      {  1,   1,   1,   1 },
      { 10,  10,  10,  10 }
    },
    {
      {  6,   3,   6,   3 },
      {  0,   0,   0,   0 },
      {  7,   7,   7,   7 },
      { 19,  19,  12,  19 }
    },
    {
      { 20,  18,  20,  18 },
      { 15,  15,  15,  15 },
      { 17,   4,  17,   4 },
      { 10,  13,  12,  13 }
    }
  },
  {
    {
      {  2,   2,  11,   2 },
      { 16,  16,  16,  16 },
      { 15,  15,  15,  15 },
      {  9,   9,  12,   9 }
    },
    {
      {  5,   8,   5,   8 },
      { 14,  14,  14,  14 },
      {  1,   1,   1,   1 },
      { 10,  10,  10,  10 }
    },
    {
      {  6,   3,   6,   3 },
      {  0,   0,   0,   0 },
      {  7,   7,   7,   7 },
      { 19,  19,  19,  19 }
    },
    {
      { 18,  18,  20,  18 },
      { 15,  15,  15,  15 },
      { 17,   4,  17,   4 },
      { 10,  13,  10,  13 }
    }
  },
  {
    {
      { 11,   2,  11,   2 },
      { 16,  16,  16,  16 },
      {  1,  15,   1,  15 },
      {  9,   9,  12,   9 }
    },
    {
      {  5,   8,   5,   8 },
      { 14,  14,  14,  14 },
      {  1,   1,   1,   1 },
      { 10,  10,  10,  10 }
    },
    {
      {  6,   3,   6,   3 },
      {  0,   0,   0,   0 },
      {  7,   7,   7,   7 },
      { 19,  19,  19,  19 }
    },
    {
      { 20,  18,   5,  18 },
      { 15,  15,  15,  15 },
      { 20,   4,  17,   4 },
      { 10,  13,  10,  13 }
    }
  },
  {
    {
      { 11,   2,  11,   2 },
      { 16,  16,  16,  16 },
      {  1,  15,   1,  15 },
      {  9,   9,  12,   9 }
    },
    {
      {  5,   8,   5,   8 },
      { 14,  14,  14,  14 },
      {  1,   1,   1,   1 },
      { 10,  10,  10,  10 }
    },
    {
      {  6,   3,   6,   3 },
      {  0,   0,   0,   0 },
      {  7,   7,   7,   7 },
      { 19,  19,  19,  19 }
    },
    {
      { 20,  18,  10,  18 },
      { 15,  15,  15,  15 },
      { 20,   4,  17,   4 },
      { 10,  13,  10,  13 }
    }
  },
  {
    {
      {  2,   2,  11,   2 },
      { 16,  16,  16,  16 },
      { 15,  15,  15,  15 },
      { 12,   9,  12,   9 }
    },
    {
      {  5,   8,   5,   8 },
      { 14,  14,  14,  14 },
      {  1,   1,   1,   1 },
      { 10,  10,  10,  10 }
    },
    {
      {  6,   3,   6,   3 },
      {  0,   0,   0,   0 },
      {  7,   7,   7,   7 },
      { 19,  19,  12,  19 }
    },
    {
      { 20,  18,  20,  18 },
      { 15,  15,  15,  15 },
      { 17,   4,  17,   4 },
      { 10,  13,  10,  13 }
    }
  },
  {
    {
      { 11,   2,  11,   2 },
      { 16,  16,  16,  16 },
      {  1,  15,   1,  15 },
      {  9,   9,  12,   9 }
    },
    {
      {  5,   8,   5,   8 },
      { 14,  14,  14,  14 },
      {  1,   1,   1,   1 },
      { 10,  10,  10,  10 }
    },
    {
      {  6,   3,   6,   3 },
      {  0,   0,   0,   0 },
      {  7,   7,   7,   7 },
      { 19,  19,  19,  19 }
    },
    {
      { 20,  18,  10,  18 },
      { 20,  15,  15,  15 },
      { 20,   4,  17,   4 },
      { 10,  13,  10,  13 }
    }
  },
  {
    {
      { 11,   2,  11,   2 },
      { 16,  16,  16,  16 },
      {  1,  15,   1,  15 },
      {  9,   9,  12,  12 }
    },
    {
      {  5,   8,   5,   8 },
      { 14,  14,  14,  14 },
      {  1,   1,   1,   1 },
      { 10,  10,  10,  10 }
    },
    {
      {  6,   3,   6,   3 },
      {  0,   0,   0,   0 },
      {  7,   7,   7,   7 },
      { 19,  19,  12,  19 }
    },
    {
      { 20,  18,  20,  18 },
      { 15,  15,  15,  15 },
      { 20,   4,  17,   4 },
      { 20,  13,  10,  13 }
    }
  }
};

static const int NUM_AMINOACID = 20;
static const int NUM_NUCLEOTIDE = 4;
static const int POS_START_CODON_DATA = 12;
static char buffer[LEN_BUFFER + 1];

static int read_one_transl_table (FILE *fp);
static int get_id_from_buffer (int *v, char *bf);
static int populate_Nuc2AATable (int* table[][4][4],
                                 int* tableStart[][4][4],
                                 char ncbieea[], char sncbieea[],
                                 char Base1[], char Base2[], char Base3[]);
static int print_Nuc2AATable (int table[][NUM_NUCLEOTIDE][NUM_NUCLEOTIDE]);

static int print_Nuc2AATable (int table[][NUM_NUCLEOTIDE][NUM_NUCLEOTIDE])
{
  int i, j, k;
  printf ("  {\n");
  for (i = 0; i < NUM_NUCLEOTIDE; i++)
    {
      printf ("    {\n");
      for (j = 0; j < NUM_NUCLEOTIDE; j++)
        {
          printf ("      {");
          for (k = 0; k < NUM_NUCLEOTIDE; k++)
            {
              printf ("%3d", table[i][j][k]);
              if (k != 3) 
                {
              printf (", ");
                }
            }
          if (j != 3) 
            {
              printf (" },\n");
            }
          else
            {
              printf (" }\n");
            }
        }
      if (i != 3) 
        {
          printf ("    },\n");
        }
      else
        {
          printf ("    }\n");
        }
    }
  printf ("  },\n");
  return EXIT_SUCCESS;
}

static int 
populate_Nuc2AATable (int* table[][4][4],
                      int* tableStart[][4][4],
                      char ncbieaa[], char sncbieaa[],
                      char Base1[], char Base2[], char Base3[])
{
  int i;
  int b1, b2, b3; 

  for (i = 0; i < NUM_CODON; i++)
    {
      b1 = INT_NUCLEOTIDE[(int) Base1[i]];
      b2 = INT_NUCLEOTIDE[(int) Base2[i]];
      b3 = INT_NUCLEOTIDE[(int) Base3[i]];
      assert (INT_AMINOACID[(int) ncbieaa[i]] >= 0);
      assert (INT_AMINOACID[(int) ncbieaa[i]] <= 20);
      table[b1][b2][b3] = INT_AMINOACID[(int) ncbieaa[i]];
      if (sncbieaa[i] == 'M') 
        {
          tableStart[b1][b2][b3] = INT_AMINOACID[(int) sncbieaa[i]];
        } 
      else 
        {
          tableStart[b1][b2][b3] = INT_AMINOACID[(int) ncbieaa[i]];
        } 
    }
  return EXIT_SUCCESS;
}

static int 
get_id_from_buffer (int *v, char *bf)
{
  char t[LEN_BUFFER];
  int r;
  r = sscanf (bf, "%s %d %s", t, v, t);
  assert (r == 3);
  return EXIT_SUCCESS;
}

static int
read_one_transl_table (FILE *fp)
{
  char c;
  int i;
  int id;
  int transl_table_id =  PSI_CODON_TABLE_1;
  char ncbieaa[NUM_CODON];
  char sncbieaa[NUM_CODON];
  char Base1[NUM_CODON];
  char Base2[NUM_CODON];
  char Base3[NUM_CODON];
  while ((c = getc (fp)) != '{' && c != '}') {};

  if (c != '{') 
    {
      return EXIT_FAILURE;
    }
  
  fgets (buffer, LEN_BUFFER, fp);
  while (strstr(buffer, " id ") == NULL)
    {
      fgets (buffer, LEN_BUFFER, fp);
    }
  get_id_from_buffer (&id, buffer);
/*  printf ("id: %d\n", id); */

  for (i = 0; i < POS_START_CODON_DATA; i++) { getc (fp); }
  for (i = 0; i < NUM_CODON; i++) { ncbieaa[i] = getc (fp); }
  fgets (buffer, LEN_BUFFER, fp);

  for (i = 0; i < POS_START_CODON_DATA; i++) { getc (fp); }
  for (i = 0; i < NUM_CODON; i++) { sncbieaa[i] = getc (fp); }
  fgets (buffer, LEN_BUFFER, fp);

  for (i = 0; i < POS_START_CODON_DATA; i++) { getc (fp); }
  for (i = 0; i < NUM_CODON; i++) { Base1[i] = getc (fp); }
  fgets (buffer, LEN_BUFFER, fp);

  for (i = 0; i < POS_START_CODON_DATA; i++) { getc (fp); }
  for (i = 0; i < NUM_CODON; i++) { Base2[i] = getc (fp); }
  fgets (buffer, LEN_BUFFER, fp);
  
  for (i = 0; i < POS_START_CODON_DATA; i++) { getc (fp); }
  for (i = 0; i < NUM_CODON; i++) { Base3[i] = getc (fp); }
  fgets (buffer, LEN_BUFFER, fp);

  /* Let's populate Nuc2AATables */
  switch (id)
    {
      case 1: transl_table_id =  PSI_CODON_TABLE_1; break;
      case 2: transl_table_id =  PSI_CODON_TABLE_2; break;
      case 3: transl_table_id =  PSI_CODON_TABLE_3; break;
      case 4: transl_table_id =  PSI_CODON_TABLE_4; break;
      case 5: transl_table_id =  PSI_CODON_TABLE_5; break;
      case 6: transl_table_id =  PSI_CODON_TABLE_6; break;
      case 9: transl_table_id =  PSI_CODON_TABLE_9; break;
      case 10: transl_table_id =  PSI_CODON_TABLE_10; break;
      case 11: transl_table_id =  PSI_CODON_TABLE_11; break;
      case 12: transl_table_id =  PSI_CODON_TABLE_12; break;
      case 13: transl_table_id =  PSI_CODON_TABLE_13; break;
      case 14: transl_table_id =  PSI_CODON_TABLE_14; break;
      case 15: transl_table_id =  PSI_CODON_TABLE_15; break;
      case 16: transl_table_id =  PSI_CODON_TABLE_16; break;
      case 21: transl_table_id =  PSI_CODON_TABLE_21; break;
      case 22: transl_table_id =  PSI_CODON_TABLE_22; break;
      case 23: transl_table_id =  PSI_CODON_TABLE_23; break;
      default: assert (0);
    }
  populate_Nuc2AATable (Nuc2AATables[transl_table_id],
                        Nuc2AATablesStart[transl_table_id],
                        ncbieaa, sncbieaa,
                        Base1, Base2, Base3);

  print_Nuc2AATable (Nuc2AATables[transl_table_id]); 
/*
  print_Nuc2AATable (Nuc2AATablesStart[transl_table_id]);
*/

  while ((c = getc (fp)) != '}') {};

  return EXIT_SUCCESS;
}

int 
num_nucleotide (int *dna, int len, int w)
{
  int i;
  int n = 0;
  int *nuc = dna;
  for (i = 0; i < len; i++) 
    {
      if (*nuc == w) 
        {
          n++;
        }
      nuc++;
    }
  return n;
}

int 
find_codon_in_dna (int *codon, int *dna, int pos)
{
  /* ACGTACGTACGTACGTACGTACGT
        p..   .p.   ..p
   */
 
  /* 'p' positioned nucleotide */
  int whichNuc = pos % 3;
  codon[whichNuc] = dna[pos];
  assert (dna[pos] == PSI_DNA_A || dna[pos] == PSI_DNA_C
          || dna[pos] == PSI_DNA_G || dna[pos] == PSI_DNA_T);

  /* next-to-'p' positioned nucleotide */
  whichNuc++;
  if (whichNuc < 3) 
    {
      pos++;
      codon[whichNuc] = dna[pos];
      assert (dna[pos] == PSI_DNA_A || dna[pos] == PSI_DNA_C
              || dna[pos] == PSI_DNA_G || dna[pos] == PSI_DNA_T);
    }
  else
    {
      whichNuc = 0;
      pos -= 2; 
      codon[whichNuc] = dna[pos];
      assert (dna[pos] == PSI_DNA_A || dna[pos] == PSI_DNA_C
              || dna[pos] == PSI_DNA_G || dna[pos] == PSI_DNA_T);
    }

  /* the next-to-'p' positioned nucleotide */
  whichNuc++;
  if (whichNuc < 3) 
    {
      pos++;
      codon[whichNuc] = dna[pos];
      assert (dna[pos] == PSI_DNA_A || dna[pos] == PSI_DNA_C
              || dna[pos] == PSI_DNA_G || dna[pos] == PSI_DNA_T);
    }
  else
    {
      whichNuc = 0;
      pos -= 2; 
      codon[whichNuc] = dna[pos];
      assert (dna[pos] == PSI_DNA_A || dna[pos] == PSI_DNA_C
              || dna[pos] == PSI_DNA_G || dna[pos] == PSI_DNA_T);
    }

  return EXIT_SUCCESS;
}

int 
read_translaton_table (const char *fn)
{
  FILE *fp = NULL;
  int c;

  fp = fopen (fn, "r");
 
  while ((c = getc (fp)) != '{') {};

  printf ("static int Nuc2AATables[17][4][4][4] = \n");
/*
  printf ("static int Nuc2AATablesStart[17][4][4][4] = \n");
*/
  printf ("{\n");
  while ((c = read_one_transl_table (fp)) == EXIT_SUCCESS) {};
  printf ("};\n");
  
  fclose (fp);
  return EXIT_SUCCESS;
}

int 
choose_transl_table (int id, int init)
{
  int i, j, k;

  init_codon = init;
  switch (id)
    {
      case 1: translation_table_id =  PSI_CODON_TABLE_1; break;
      case 2: translation_table_id =  PSI_CODON_TABLE_2; break;
      case 3: translation_table_id =  PSI_CODON_TABLE_3; break;
      case 4: translation_table_id =  PSI_CODON_TABLE_4; break;
      case 5: translation_table_id =  PSI_CODON_TABLE_5; break;
      case 6: translation_table_id =  PSI_CODON_TABLE_6; break;
      case 9: translation_table_id =  PSI_CODON_TABLE_9; break;
      case 10: translation_table_id =  PSI_CODON_TABLE_10; break;
      case 11: translation_table_id =  PSI_CODON_TABLE_11; break;
      case 12: translation_table_id =  PSI_CODON_TABLE_12; break;
      case 13: translation_table_id =  PSI_CODON_TABLE_13; break;
      case 14: translation_table_id =  PSI_CODON_TABLE_14; break;
      case 15: translation_table_id =  PSI_CODON_TABLE_15; break;
      case 16: translation_table_id =  PSI_CODON_TABLE_16; break;
      case 21: translation_table_id =  PSI_CODON_TABLE_21; break;
      case 22: translation_table_id =  PSI_CODON_TABLE_22; break;
      case 23: translation_table_id =  PSI_CODON_TABLE_23; break;
      default: assert (0);
    }

  for (i = 0; i < NUM_NUCLEOTIDE; i++)
    {
      for (j = 0; j < NUM_NUCLEOTIDE; j++)
        {
          for (k = 0; k < NUM_NUCLEOTIDE; k++)
            {
              Nuc2AATable[i][j][k] = 
                Nuc2AATables[translation_table_id][i][j][k];
              Nuc2AATableStart[i][j][k] = 
                Nuc2AATablesStart[translation_table_id][i][j][k];
            }
        }
    }
  return EXIT_SUCCESS;
}

int 
dna_char2int (char *str_dna, int *dna, int len)
{
  int i;
  int j;
  char nucleotide;
  for (i = 0; i < len; i++) 
    { 
      nucleotide = str_dna[i];
      assert (nucleotide == 'A' 
              || nucleotide == 'C'
              || nucleotide == 'G'
              || nucleotide == 'T');
      j = (int) nucleotide;
      dna[i] = INT_NUCLEOTIDE[j];
    }
  return EXIT_SUCCESS;
}

int dna_int2char (int *dna, char *str_dna, int len)
{
   int i;
   int j;
   for (i = 0; i < len; i++) {
      j = dna[i];
      assert (j == PSI_DNA_A
              || j == PSI_DNA_C
              || j == PSI_DNA_G
              || j == PSI_DNA_T);
      str_dna[i] = NUCLEOTIDE[j];
   }
   str_dna[len] = '\0';
   return EXIT_SUCCESS;
}

int 
pro_char2int (char *str_pro, int *pro, int len)
{
   int i;
   int j;
   char aa;
   for (i = 0; i < len; i++) {
      aa = str_pro[i];
      /* TODO: assert if aa is not a valid amino acid letter */
      j = (int) aa;
      pro[i] = INT_AMINOACID[j];
   }
   return EXIT_SUCCESS;
}

int 
pro_int2char (int *pro, char *str_pro, int len)
{
   int i;
   int j;
   for (i = 0; i < len; i++) {
      j = pro[i];
/*
fprintf (stderr, "%d", j);
*/
      if (j >= 0 && j <= 20) {
         str_pro[i] = AMINOACID[j];
      } else {
         return 1; 
      }
   }
   str_pro[len] = '\0';
/*
fprintf (stderr, "\n");
*/
   return EXIT_SUCCESS;
}

int 
dna2protein (int *dna, int *protein, int len)
{
  int i;
  int aa;

  if (init_codon == 1)
    protein[0] = Nuc2AATableStart[dna[0]][dna[1]][dna[2]];
  else
    protein[0] = Nuc2AATable[dna[0]][dna[1]][dna[2]];

/*
  if (protein[0] == PSI_AA_STP)
    return EXIT_FAILURE;
*/

  for (i = 1; i < len; i++) 
    {
      aa = Nuc2AATable[dna[3*i]][dna[3*i + 1]][dna[3*i + 2]];
/*
      if (aa == PSI_AA_STP)
        return EXIT_FAILURE;
*/
      protein[i] = aa;
    }

  for (i = 0; i < len; i++) 
    {
      if (protein[i] == PSI_AA_STP)
        return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}

/*
 * pi: double array of size four - already allocated
 * c : double array of size 20 - already allocated
 */
int
codon_frequency (double *pi, double *c) 
{
   int i, j, k;
   if (pi == NULL || c == NULL) {
      return 1; /* ERR_MEMORY; */
   }
   for (i = 0; i < NUM_AMINOACID; i++) {
      c[i] = 0;
   }
   for (i = 0; i < NUM_NUCLEOTIDE; i++) {
      for (j = 0; j < NUM_NUCLEOTIDE; j++) {
         for (k = 0; k < NUM_NUCLEOTIDE; k++) {
            if (Nuc2AATable[i][j][k] != 20) {
               /* c[Nuc2AATable[i][j][k]] += logsum (log(pi[i]*pi[j]*pi[k]),
                                                  c[Nuc2AATable[i][j][k]]); */
               c[Nuc2AATable[i][j][k]] += (pi[i]*pi[j]*pi[k]);
            }
         }
      }
   }

   for (i = 0; i < NUM_AMINOACID; i++) {
      c[i] = log(c[i]);
   }

   return EXIT_SUCCESS; 
}

int
aa_codon_is (int codon[])
{
  return Nuc2AATable[codon[0]][codon[1]][codon[2]];
}

int
change_translation_table(int i)
{
  psi_fatal ("no implementaton of change_translation_table in seq.c/h");
  int r = EXIT_SUCCESS;
  switch (i) {
  case 1:
     /* No Code */
     break;
  case 4:
     Nuc2AATable[PSI_DNA_T][PSI_DNA_G][PSI_DNA_A] = INT_AMINOACID[(int)'W'];
     break;
  case 6:
     Nuc2AATable[PSI_DNA_T][PSI_DNA_A][PSI_DNA_A] = INT_AMINOACID[(int)'Q'];
     Nuc2AATable[PSI_DNA_T][PSI_DNA_A][PSI_DNA_G] = INT_AMINOACID[(int)'Q'];
     break;
  default:
     r = 1; /* ERR_TRANSLATION; */
  }
  return r;
}

int
psi_seq_diff_seqs (int *seq_a, int *seq_b, int len)
{
  int i, n;
  n = 0;
  for (i = 0; i < len; i++)
    {
      if (seq_a[i] != seq_b[i])
        n++;
    }
  return n;
}

double seq_put_out (double d)
{
  return d;
}
