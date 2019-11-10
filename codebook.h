// mapping of ranges of harmonic energies h1-h15 to phonetic interpretation
// each entry is the acceptable range from minimum and maximum inclusive
// Energy range is 0 to 15, represeted in hex format, 0x[min][max]

#define CBSIZE 30
const unsigned char cb[CBSIZE][16]={

// silence
0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,'_',

// 'ah' sound in 'i' (ae), 'j' (jae), 'r' (ar), 'y' (wae).
0x7c,0x7c,0x8f,0x8e,0x6b,0x48,0x25,0x25,0x37,0x4a,0x48,0x48,0x49,0x4a,0x48,'a',

// 'ay' sound used in 'a' (A), 'h' (AC), 'k' (kA)
0x3e,0x3e,0x26,0x04,0x04,0x05,0x17,0x4b,0x8e,0x8e,0x7d,0x7b,0x6c,0x6c,0x5f,'A',

// 'b' sound used in 'b' (b), 'p' (bkE) and 'w' (dublyO)
0x3d,0x39,0x26,0x26,0x39,0x5d,0x6e,0x5d,0x6c,0x5e,0x5c,0x3a,0x3a,0x49,0x4a,'b',

// 'd' used in 'd' (dE), 'w' (dublO), 'z' (sed)
0x07,0x29,0x39,0x39,0x29,0x3a,0x4a,0x4b,0x5c,0x6c,0x7c,0x6d,0x6c,0x6e,0x6e,'d',

// 'ee'and 'y' sound used in 'b' (bE), 'c' (sE), 'd' (dE), 'e' (E), 'g' (jE),
//                    'p' (bkE), 't' (tE), 'v' (vE or fE), 'u' (EO)
0xaf,0x49,0x15,0x04,0x04,0x06,0x16,0x3a,0x5b,0x7d,0x6d,0x7d,0x7d,0x6e,0x6d,'E',

// 'eh' sound used in 'f' (ef), 'l' (el), 'm' (em), 'n' (en), 
//                   's' (es), 'x' (eks), 'z' (sed)
0x4c,0x5b,0x4d,0x26,0x26,0x3a,0x6c,0x5b,0x5a,0x6c,0x49,0x47,0x59,0x5a,0x49,'e',
// 'f' sound from 'f' (ef),  'v' (fE).
0x15,0x24,0x14,0x35,0x56,0x69,0x6b,0x69,0x6a,0x8d,0x8b,0x79,0x6a,0x7a,0x7b,'f',
0xef,0x9b,0x89,0x77,0x67,0x69,0x6a,0x69,0x58,0x56,0x66,0x55,0x45,0x45,0x45,'f',
0x99,0x99,0x99,0x66,0x55,0x69,0xbb,0x89,0x78,0xcc,0x66,0x44,0x55,0x88,0x44,'f',
0x11,0x22,0x33,0x77,0x88,0x88,0x99,0x88,0x77,0x99,0xaa,0x88,0x99,0xbb,0x88,'f',

// 'j/ch' sound from 'g' (jE), 'h' (Aj), 'j' (jae).
0x01,0x01,0x02,0x13,0x14,0x25,0x49,0x5b,0x8d,0xaf,0xcf,0xbf,0xaf,0xaf,0xae,'j',
0x33,0x44,0x55,0x55,0x66,0x77,0x77,0x88,0x88,0x99,0x99,0xaa,0xaa,0x99,0x88,'j',

// 'k' sound used in 'k' (kA), 'p' (bkE), 'q' (KyO).
0x05,0x27,0x25,0x14,0x14,0x26,0x47,0x5a,0x6b,0x6d,0xaf,0x9e,0x8e,0x8e,0x9d,'k',

// 'l' sound in 'l' (el)
0x8d,0x7d,0x7d,0x6d,0x48,0x36,0x26,0x27,0x27,0x47,0x48,0x59,0x59,0x59,0x5c,'l',

// 'm' sound in 'm' (em)
0xac,0x56,0x35,0x36,0x59,0x46,0x45,0x58,0x7a,0x79,0x79,0x69,0x68,0x78,0x69,'m',

// 'n' sound in 'n' (en)
0xac,0x57,0x35,0x25,0x36,0x38,0x47,0x5b,0x78,0x7a,0x7e,0x68,0x58,0x69,0x7a,'n',

// 'o' sound in 'o' (ow),
0xcf,0xbf,0xbf,0x79,0x45,0x34,0x34,0x35,0x35,0x46,0x57,0x58,0x48,0x59,0x56,'o',

// 'oo' sound in 'q' (kO), 'u' (EO), 'w' (dublyO)
0xae,0x5d,0x36,0x24,0x14,0x26,0x6a,0x6d,0x69,0x8f,0x5a,0x4a,0x59,0x5b,0x48,'O',

// 'r' sound only used in 'r' (ar).
0x7b,0x6b,0x6c,0x8c,0x7d,0x7b,0x7b,0x5a,0x37,0x25,0x17,0x27,0x47,0x37,0x47,'r',

// 's' sound in 'c' (sE), 's' (es), 'x' (e_s), 'z' (sed)
0x03,0x03,0x05,0x14,0x25,0x37,0x5b,0x59,0x6b,0x8c,0x9e,0x9f,0xaf,0xcf,0xcf,'s',

// plosive 't' sound in 't' (tE)
0x01,0x02,0x14,0x25,0x37,0x3a,0x59,0x58,0x7c,0x8c,0xac,0xbf,0x9e,0xad,0xac,'t',
0x01,0x22,0x44,0x66,0x77,0x77,0x66,0x77,0x88,0x99,0x99,0xbb,0xbb,0x99,0xaa,'t',

// 'u' sound only present in 'w' (dublyO)
0x9e,0xaf,0x7e,0x6a,0x6a,0x5a,0x3a,0x27,0x16,0x27,0x58,0x58,0x48,0x57,0x38,'u',

// 'v' sound in 'v' (vE or fE)
0x9f,0x48,0x36,0x24,0x45,0x58,0x7a,0x68,0x68,0x68,0x5c,0x69,0x68,0x69,0x6a,'v',
0x58,0x34,0x23,0x24,0x26,0x39,0x58,0x69,0x79,0x78,0x8b,0x8a,0x8a,0x9b,0x8c,'v',
0x03,0x13,0x26,0x36,0x35,0x47,0x7b,0x7c,0x7b,0x8a,0x8d,0x9b,0x8b,0x8c,0x8c,'v',

// 'w' sound in 'y' (wae)
0xaf,0x8f,0x7f,0x6e,0x49,0x36,0x27,0x27,0x27,0x37,0x38,0x38,0x39,0x49,0x38,'w',
0x66,0x88,0xab,0x9b,0x68,0x56,0x56,0x77,0x67,0x46,0x57,0x77,0x68,0x68,0x68,'w',
0x99,0x99,0xcc,0xbb,0x88,0x66,0x44,0x44,0x55,0x66,0x66,0x66,0x66,0x66,0x77,'w',

};