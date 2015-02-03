/*   mmdblocl.c
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*            National Center for Biotechnology Information (NCBI)
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government do not place any restriction on its use or reproduction.
*  We would, however, appreciate having the NCBI and the author cited in
*  any work or product based on this material
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
* ===========================================================================
*
* File Name:  mmdblocl.c
*
* Author:  Christopher Hogue
*
* Version Creation Date:  14 January 97  
*
* $Revision: 1.9 $
*
* File Description: Used to provide Biostruc data to MMDB-API when
*  MMDB files are on a local filesystem in standard MMDB ftp-site format.
*                   
* Modifications:  
* --------------------------------------------------------------------------
* $Log: mmdblocl.c,v $
* Revision 1.9  2001/09/10 16:49:38  michel
* added SLRI_MMDBBiostrucGet to call MMDBGetBiostrucGet (for clustal library)
*
* Revision 1.8  2001/06/21 16:19:10  feldman
* Fixed potential fatal bug
*
* Revision 1.7  2001/06/19 20:34:38  feldman
* Allow gzip to work on WINDOWS as well
*
* Revision 1.6  2001/03/29 20:56:59  feldman
* Fixed some Windows compiler warnings, pertaining to typecasting, mostly
* doubles to floats
*
* Revision 1.5  2001/01/12 20:01:51  feldman
* Added functionality to maketrj to call RPSBlast and then look up
* results in local CDD database - this functionality makes aligntraj
* obsolete, so that program will not be developed further
*
* Major functions were cut and paste from cddserver.c in toolkit and
* adapted into clust.c.
* Note this addition required 2 more parameters in the .foldtrajrc file
* for RPSBlast database name and CDD path
*
* Revision 1.4  2000/10/06 17:34:25  feldman
* Fixed potential problem on SGI (again)
*
* Revision 1.3  2000/10/03 16:59:50  feldman
* reverted to original file, must find an alternate solution for SGI version..
*
* Revision 6.0  1997/08/25 18:11:25  madden
* Revision changed to 6.0
*
* Revision 1.5  1997/07/14 18:29:03  madej
* Changed gunzip access to get pathname from config file.
*
 *
 * Revision 1.4  1997/03/25  22:24:42  hogue
 * Fixed an off-by-one error.
 *
 * Revision 1.3  1997/03/19  16:56:52  hogue
 * Added MMDB Iterators, made better error messages
 *
 * Revision 1.2  1997/02/05  23:13:56  hogue
 * Added gunzip as default method of decompression of files on UNIX systems.
 * Directs gunzip to a temp file, then deletes it after the load.
 *
 * Revision 1.1  1997/01/27  18:42:02  hogue
 * Initial revision
 *
*
*
* ==========================================================================
*/

/*********************************************************************************/
/* This file abstracts the calls to Network Entrez that were previously          */
/* embedded in mmdbapi1.c  Linking with this file forces local MMDB storage      */
/* Local MMDB storage is defined as that distributed on the MMDB FTP site at     */
/* ftp://ncbi.nlm.nih.gov/pub/mmdb/mmdb/Biostruc                                 */
/* which contains MMDB ASN.1 files of the form 12345.val and an index "mmdb.idx" */
/* This is set up so that one can mirror the MMDB database from the ftp site.    */

/********   NOTE - .mmdbrc (UNIX) mmdb.ini (Win) mmdb.cnf (Mac) required   *******/
/* containing paths to database and index file name:                         *****/
/* 
[MMDB]
Database = /net/dorothy/export/home/hogue/mmdb/Oct96/database/
Index    = mmdb.idx
 */



#include <ncbi.h>
#include <mmdbapi1.h>
#include <mmdbdata.h>
#include <mmdblocl.h>
#include <accentr.h>
#include <accutils.h>

/*#ifdef OS_UNIX*/
/* works on windows too, with gzip.exe installed */
#define MMDB_UNIXCOMPRESSED 1
/*#endif*/

/* LoadMMDBIdx reads the local MMDB index "mmdb.idx" into memory */
/* index format is
3\n
1 1ABC\n
3 3INS\n
234567 4GWA\n
eof
*/


static Int4Ptr pI4Mmdbidx = NULL;
static Uint4Ptr pU4Pdbidx = NULL;
static long iMMDBSize = 0;
static Char database[PATH_MAX];
static Char indexname[256];
static Int4Ptr pI4position = NULL;
static Char gunzip[PATH_MAX];
static Char pdbCode[5];

static Uint4 HashPDBCode(CharPtr pcPdb)
{

  Uint4 iPdbhash;
  Uint4 c1, c2, c3, c4;

  if (!pcPdb) return 0;
  if (StringLen(pcPdb) != 4) return 0;
  c1 =  pcPdb[3];
  c2 =  pcPdb[2];
  c3 =  pcPdb[1];
  c4 =  pcPdb[0];
  iPdbhash =  c1 + (c2 << 8) + (c3 << 16) + (c4 << 24);
  
  return iPdbhash;
}

 
Int4 LIBCALL GetMMDBSize(void)
{
  return (Int4) iMMDBSize;  
}

void LIBCALL ResetMMDBIterate(void)
{
  pI4position = pI4Mmdbidx;
}

void LIBCALL SetMMDBIterateStart(Int4 MMDBid)
{
  Int4 i;
  Int4Ptr      pI4Mmdb;
  
  pI4Mmdb = pI4Mmdbidx;
  for (i=0; i < (iMMDBSize) ; i++)
      {
        if (MMDBid == *pI4Mmdb)
	 {
	      pI4position = pI4Mmdb;
              return;
	 }
	pI4Mmdb++;
      }  
  pI4position = pI4Mmdbidx;;
  return;
}

Int4 LIBCALL GetMMDBIdCurrent(void)
{
  return *pI4position;
}


Boolean LIBCALL IterateMMDBId(void)
{
  
  pI4position++;
  if (pI4position == (Int4Ptr) &pI4Mmdbidx[iMMDBSize])
   {
     pI4position = pI4Mmdbidx;
     return FALSE;
   }
  return TRUE;
}


static CharPtr GetHashedPDBCode(Uint4 ulPDB)
{
 
  Uint4 uitemp1, uitemp2, uitemp3, uitemp4;
  
  uitemp1 = ulPDB & (0xFF000000);
  uitemp2 = ulPDB & (0x00FF0000);
  uitemp3 = ulPDB & (0x0000FF00);
  uitemp4 = ulPDB & (0x000000FF);
  
  pdbCode[0] = (Char) (uitemp1 >> 24);
  pdbCode[1] = (Char) (uitemp2 >> 16);
  pdbCode[2] = (Char) (uitemp3 >> 8);
  pdbCode[3] = (Char) (uitemp4);  
  pdbCode[4] = '\0';

  return (CharPtr) pdbCode;
}


static Boolean LoadMMDBIdx(CharPtr db,  CharPtr idx)
{

   FILE *f;
   Char fullpath [PATH_MAX];
   Char pcBuf[100];
   CharPtr pcTest;
   Char  pcMmdb[20];
   Char  pcPdb[20];
   Uint4Ptr     pU4Pdb;
   Int4Ptr      pI4Mmdb;
   long count = 0;
   long mmdbid = 0;
/* StrToLong stuff */
  Nlm_Char     ch;
  Nlm_Int2     i;
  Nlm_Int2     len;
  Nlm_Char     local [64];
  Nlm_Boolean  nodigits;
  Nlm_Boolean  rsult;
  long int     val;
 

   StringCpy(fullpath, db);
   StringCat(fullpath, idx);
   if ((f = FileOpen (fullpath, "r")) == NULL) 
  	 {
	    	ErrPostEx(SEV_FATAL,0,0, "MMDBInit - Failed to load MMDB index - check config file.");
		return FALSE;
	  } 
   
   pcTest = fgets(pcBuf, (size_t)100, f);
   if (pcTest)
     {
       sscanf(pcBuf, "%ld", &iMMDBSize);
     }
   if ((iMMDBSize == 0) || (iMMDBSize > 100000))
     {
	    	ErrPostEx(SEV_FATAL,0,0, "Internal - LoadMMDBIdx() Failure 2;");
		return FALSE;
     } 
   pI4Mmdbidx = (Int4Ptr) MemNew((size_t) ((iMMDBSize)*sizeof(Int4)));
   pU4Pdbidx = (Uint4Ptr) MemNew((size_t) ((iMMDBSize)*sizeof(Uint4)));

   if ((!pI4Mmdbidx) || (!pU4Pdbidx))
     {
	    	ErrPostEx(SEV_FATAL,0,0, "Internal - LoadMMDBIdx() Out of Memory;");
		return FALSE;
     }
   
   pI4Mmdb = pI4Mmdbidx;
   pU4Pdb = pU4Pdbidx;
   
   do 
     {
      pcBuf[0] = '\0';
      pcTest = fgets(pcBuf,  (size_t)100,  f);
      if (pcTest)
        {
	    sscanf(pcBuf, "%s%s",  pcMmdb ,  pcPdb);
	    StrUpper(pcPdb);	 
	    mmdbid = 0;
 	    
/* this is code from Nlm StrToLong */
  rsult = FALSE;
  len = (Nlm_Int2) Nlm_StringLen (pcMmdb);
  if (len != 0) {
    rsult = TRUE;
    nodigits = TRUE;
    for (i = 0; i < len; i++) {
      ch = pcMmdb [i];
      if (ch == ' ' || ch == '+' || ch == '-') {
      } else if (ch < '0' || ch > '9') {
        rsult = FALSE;
      } else {
        nodigits = FALSE;
      }
    }
    if (nodigits) {
      rsult = FALSE;
    }
    if (rsult) {
      Nlm_StringNCpy (local, pcMmdb, sizeof (local));
      if (sscanf (local, "%ld", &val) == 1) {
        mmdbid = val;
      }
    }
  }


/* printf("%ld %s %ld %ld\n", mmdbid, pcPdb, (long) count, (long) iMMDBSize-1);  */
            if ((mmdbid == 0) ||  (StringLen(pcPdb) > 4))
	      {
                MMDBFini();
	    	ErrPostEx(SEV_FATAL,0,0, "Internal - LoadMMDBIdx() Bad Index Values");
		return FALSE;
	       }
    
            *pI4Mmdb = mmdbid;
	    *pU4Pdb =  HashPDBCode(pcPdb);
	    pI4Mmdb++;
	    pU4Pdb++;
	    count++;
	}   
       if (count > (iMMDBSize))
       {
                MMDBFini();
	    	ErrPostEx(SEV_FATAL,0,0, "Internal - LoadMMDBIdx() Index count is wrong at top of file;");
		return FALSE;
       }
     } while (pcTest);   
   FileClose(f);
   return TRUE;
}



static void RemoveMMDBIdx()
{
 
   if (pI4Mmdbidx) pI4Mmdbidx=MemFree(pI4Mmdbidx);
   if (pU4Pdbidx) pU4Pdbidx=MemFree(pU4Pdbidx);

  return;
}




CharPtr  LIBCALL PDBidFromMMDB (Int4 iMMDB)
{
  register Uint4Ptr pU4Pdb = NULL;
  register Int4Ptr pI4Mmdb = NULL;
  Int4 i;
  
/* uses the index to make a lookup */
 
  ASSERT ((pI4Mmdbidx != NULL) && (pU4Pdbidx != NULL));
  
  pI4Mmdb = pI4Mmdbidx;
  pU4Pdb = pU4Pdbidx;
  for (i=0; i < iMMDBSize ; i++)
    {
        if (iMMDB == *pI4Mmdb)
	 {
	       return GetHashedPDBCode(*pU4Pdb);
	 }
	pU4Pdb++;
	pI4Mmdb++;
    } 
 return NULL;
}


CharPtr LIBCALL GetPDBIdCurrent(void)
{
    return PDBidFromMMDB (*pI4position);
}


DocUid  LIBCALL MMDBEvalPDB (CharPtr str)
{
  register Uint4Ptr pU4Pdb = NULL;
  register Int4Ptr pI4Mmdb = NULL;
  register Uint4   iHash = 0;
  Char eval[10];
  long i;
  
/* uses the index to make a lookup */
 

  ASSERT ((pI4Mmdbidx != NULL) && (pU4Pdbidx != NULL));
  
  pI4Mmdb = pI4Mmdbidx;
  pU4Pdb = pU4Pdbidx;
  StringNCpy(eval,  str,  4);
  eval[4] = '\0';
  StrUpper(eval);
  iHash = HashPDBCode(eval);
  for (i=0; i < iMMDBSize ; i++)
    {
        if (iHash == *pU4Pdb)
	 {
	      return (DocUid) *pI4Mmdb;
	 }
	
	pU4Pdb++;
	pI4Mmdb++;
    } 
 return (DocUid) 0;
}

 


BiostrucPtr LIBCALL MMDBBiostrucGet (DocUid uid, Int4 mdlLvl, Int4 maxModels)
{
   AsnIoPtr aip = NULL;
   Char path[PATH_MAX];
#ifdef MMDB_UNIXCOMPRESSED
   Char compath[PATH_MAX];
   Char tempfile[PATH_MAX];
#endif
   Char pcId[20];    
   Int4 iFileExists = 0;
   int  iAvail = 1;
   BiostrucPtr pbs = NULL;

   sprintf(pcId, "%ld", (long) uid);
   path[0] = '\0';
   StringCpy(path, database);
   StringCat(path, pcId);
   StringCat(path, ".val");

#ifdef MMDB_UNIXCOMPRESSED
   compath[0] = '\0';
   StringCpy(tempfile, TmpNam(NULL)); 
   sprintf(compath,"%s -c %s.gz > %s\n", gunzip, path, tempfile);
   iAvail = system(compath);
   if (iAvail != 0)
      {
 	ErrPostEx(SEV_FATAL,0,0, "MMDBBiostrucGet failed: Can't find gunzip in path.\n");
	return NULL;
      }
   iFileExists = FileLength(tempfile);
   if (iFileExists == 0)
      {
        return NULL;
      }
      
    aip = AsnIoOpen(tempfile, "rb");
#else
   iFileExists = FileLength(path);
   if (iFileExists == 0)
      {
        return NULL;
      }
      
    aip = AsnIoOpen(path, "rb");
#endif

  
    if (aip) 
     {
       pbs = BiostrucAsnGet(aip, NULL,  mdlLvl,  maxModels);
/*       pbs = BiostrucAsnRead(aip, NULL);*/
       AsnIoClose (aip);
     }

#ifdef MMDB_UNIXCOMPRESSED 
    FileRemove(tempfile);
#endif

    if (!pbs) return NULL;  
    return pbs;
} 


BiostrucPtr LIBCALL SLRI_MMDBBiostrucGet (DocUid uid, Int4 mdlLvl, Int4 maxModels)
{
	return MMDBBiostrucGet (uid, mdlLvl, maxModels);
}




Boolean LIBCALL MMDBInit (void)
{

 

  GetAppParam("mmdb", "MMDB", "Database", "", database, PATH_MAX);
  if (database[0] == 0)
    {
      	ErrPostEx(SEV_FATAL,0,0, "MMDBInit failed: config file / path to data missing.\n");
	return FALSE;
    }
  GetAppParam("mmdb", "MMDB", "Index", "", indexname, (size_t) 256 * (sizeof(char)));
  if (indexname[0] == 0)
    {
      	ErrPostEx(SEV_FATAL,0,0, "MMDBInit failed: config file / index file missing.\n");
	return FALSE;
    }
  GetAppParam("mmdb", "MMDB", "Gunzip", "", gunzip, (size_t) 256 * (sizeof(char)));
  if (gunzip[0] == 0)
    {
      	ErrPostEx(SEV_FATAL,0,0, "MMDBInit failed: gunzip path missing.\n");
	return FALSE;
    }

  return LoadMMDBIdx(database,  indexname);

}


void LIBCALL MMDBFini (void)
{

  RemoveMMDBIdx();
  return;
}
