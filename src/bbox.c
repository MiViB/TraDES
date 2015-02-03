/* 
Portions Copyright (c) 2007-2012 Hogue Laboratory, National University of Singapore
Portions Copyright (c) 1997-2007 Hogue Laboratory, University of Toronto
Portions Copyright (c) 1997-2005 Hogue Laboratory, Mount Sinai Hospital, Toronto


All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met: 

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer. 

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution. 


THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

THE COPYRIGHT HOLDERS HAVE NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, 
UPDATES, ENHANCEMENTS OR MODIFICATIONS.

Principal Contact: Christopher W.V. Hogue  chogue@blueprint.org

Contributing Hogue Laboratory Developers and Institutional Collaborators: 
(SLRI and TraDES Directories)
Christopher Hogue (mmdbapi, b-d tree, Cn3D, trades, solvate, valmerge, str2trj, seq2trj). 
Howard Feldman (salign, base foldtraj & maketraj libraries)
Contributions from John J. Salama (VisTraj) and Kevin Snyder (valmerge)
Phillipe Phan (analyzeMovie) Adrian Heilbut, Mark Kotowycz, Van Le (trj visualization),
Michel Dumontier (scoring function libraries, statistical functions) 
Elena Garderman (software maintenance), Mingxi Yao (valtopdb),
Gil Alteroviz (GOR implementation),  Boris Steipe (University of Toronto) (fragment based construction)
Brendan McConkey (University of Waterloo) & Michael Brougham(VSCORE scoring function)
Jonathan Kans (NCBI) (vistraj)

*/


/* routines for finding bounding boxes around points and doing
   certain calculations with them */

/* includes */
#include <mmdbtraj.h>

/* arbitrary value to indicate BD crash */
#define EQUAL -1	

/* global variables */
static ValNodePtr WorldsHead=NULL;
static ValNodePtr vnpNearAtoms;
static ValNodePtr vnpNearAtomsLast;
static vec minbound;
static vec maxbound;
static FloatLo circradsq;
static vec sumbound;

/* takes any PMSD, PMMD or PMGD (cast to a PFB) in a given Model, and
   finds the smallest bounding box containing that object, returning 1
   upon success or 0 if some error occurs; result must be an array
   declared by the user like so:

   vec bigbox[2];
   ....
   error=BoundBox(pfb,modelnum,bigbox)

   Upon returning, bigbox[0] will be a vector with co-ordinates at one
   corner of the box (minumum x, y and z) while bigbox[1] will be at the
   opposite corner (maximum x, y and z) */
Int2 BoundBox(PFB pfbThis, Int2 Model, vec *result)
{
	vec min,max;
	PMAD pmadHere;
	PVNAL pvnalHere;
	PALD paldHere;	
	FloatLoPtr Data;
	ValNodePtr vnpHere,vnpHead;
	Byte cnt,Set=0;
	Uint1 AtomicNo;
	FloatLo AtomicRadius;

	if (pfbThis==NULL) return 0;
	/* can only use on MSD, MMD or MGD node */
	if (!(pfbThis->bMe & (AM_MSD | AM_MMD | AM_MGD))) return 0;
	/* use ConvertNode to get a list of atoms in the object */
	vnpHead=ConvertNode(pfbThis,AM_MAD);
	vnpHere=vnpHead;
	/* walk down atoms list */
	while(vnpHere) {
		pmadHere=(PMAD)(vnpHere->data.ptrvalue);
		if (pmadHere!=NULL) {
			/* determine atomic # and hence radius where
			   possible, otherwise use 0 */
			/* see bbox.h for definition of "radius" */ 
			AtomicNo=(pmadHere->pvnmaLink)->choice;
			AtomicRadius=ElementSize(AtomicNo);
			/* walk through atom locations */
			pvnalHere=pmadHere->pvnalLocate;
			while (pvnalHere) {
				/* until desired model found */
				if (pvnalHere->choice==Model) {
					/* handle degenerate co-ords */
					paldHere=(PALD)(pvnalHere->data.ptrvalue);
					while(paldHere) {
						Data=paldHere->pflvData;
						/* if first atom, set it
						   as the bound */
						if (!Set) {
							Set=1;
							for(cnt=0;cnt<3;cnt++) {
								min[cnt]=Data[cnt]-AtomicRadius;
								max[cnt]=Data[cnt]+AtomicRadius;
							}
						}
						/* else compare to largest
						   so far, taking atomic
						   radius into account */ 
						else for(cnt=0;cnt<3;cnt++) {
							if (Data[cnt]+AtomicRadius>max[cnt])
								max[cnt]=Data[cnt]+AtomicRadius;
							if (Data[cnt]-AtomicRadius<min[cnt])
								min[cnt]=Data[cnt]-AtomicRadius;
						}
						paldHere=paldHere->next;
					}
				}
				pvnalHere=pvnalHere->next;
			}
		}
		vnpHere=vnpHere->next;
	}
	/* free list generated by ConvertNode */
	ValNodeFree(vnpHead);
	/* save result */
	result[0][0]=min[0];	
	result[0][1]=min[1];	
	result[0][2]=min[2];	
	result[1][0]=max[0];	
	result[1][1]=max[1];	
	result[1][2]=max[2];
	return 1;
}

/* acts indentically to BoundBox above except takes a world as its
   argument */
Int2 BoundBoxWorld(PWS pwsWorld, vec *finalbox)
{
	vec bigbox[2];
	Int2 error=1,flag=0,cnt;

	/* walk through each entry in the world */
	while(pwsWorld) {
		if (pwsWorld->pfbThis!=NULL) {
			error=BoundBox(pwsWorld->pfbThis,pwsWorld->Model,bigbox);
			if (!error) return 0;
			/* compare to largest bound so far and
			   update, keeping final box in finalbox */
			for (cnt=0;cnt<3;cnt++) {
				if ((!flag) || (bigbox[0][cnt]<finalbox[0][cnt]))
					finalbox[0][cnt]=bigbox[0][cnt];
				if ((!flag) || (bigbox[1][cnt]>finalbox[1][cnt]))
					finalbox[1][cnt]=bigbox[1][cnt];
			}
			flag=1;
		}
		pwsWorld=pwsWorld->next;
	}
	/* success */
	return 1;
}

/* takes any PMSD, PMMD or PMGD (cast as a PFB) in a given Model and
   returns the radius of a sphere enclosing that object; centre
   must be a vec declared by the user, and will contain a vector with
   co-ordinates at the centre of this sphere upon returning from the
   function; note: this may not be the smallest possible sphere, but
   will usually be close - the radius will be at most sqrt(3) larger
   than it needs to be (when furthest point lies at the centre of a
   face of the bounding cube) */
FloatLo BoundSphere(PFB pfbThis, Int2 Model, vec centre)
{
	Int2 err,cnt;
	vec bigbox[2];
	FloatLo radius,radtemp;
	
	/* first get the bounding box */
	err=BoundBox(pfbThis,Model,bigbox);
	if(!err)
		return 0;
	radius=0;
	for (cnt=0;cnt<3;cnt++) {
		/* centre of sphere is just centre of box */
		centre[cnt]=(bigbox[0][cnt]+bigbox[1][cnt])*0.5;
		radtemp=bigbox[1][cnt]-centre[cnt];
		/* radius of sphere enclosing box=sqrt(a^2+b^2+c^2) where
		   box has dimensions 2a x 2b x 2c (Pythagoras) */
		radius=radius+radtemp*radtemp;
	}
	return sqrt(radius);
}

/* acts identically to BoundSphere above, but takes a world as its
   argument */
FloatLo BoundSphereWorld(PWS pwsWorld, vec centre)
{
	Int2 err,cnt;
	vec bigbox[2];
	FloatLo radius,radtemp;
	
	/* first get the bounding box */
	err=BoundBoxWorld(pwsWorld,bigbox);
	if(!err)
		return 0;
	radius=0;
	for (cnt=0;cnt<3;cnt++) {
		/* centre of sphere is just centre of box */
		centre[cnt]=(bigbox[0][cnt]+bigbox[1][cnt])*0.5;
		radtemp=bigbox[1][cnt]-centre[cnt];
		/* radius of sphere enclosing box=sqrt(a^2+b^2+c^2) where
		   box has dimensions 2a x 2b x 2c (Pythagoras) */
		radius=radius+radtemp*radtemp;
	}
	return sqrt(radius);
}

/* return true if node already found in tree (not thoroughly tested..) */
Boolean BDLookup(PBDNode newnode,PBDNode treeplace,Int2 dimension)
{
	Int2 branch=0;
	Int2 d;
	FloatLo diff;

	if (newnode==NULL) return FALSE;
	for(d=dimension;d<=2;d++) {
		diff=newnode->dim[d]-treeplace->dim[d];
		if (diff!=0.0) {
			branch=(d+d+(diff>0));
			break;
		}
	}
	if (d==3)
		return TRUE;
	dimension=branch>>1;
	if(treeplace->way[branch]==NULL) {
		return FALSE;
	}
	else {
		treeplace=treeplace->way[branch];
		return BDLookup(newnode,treeplace,dimension);
	}
}

/* used internally to compare two nodes of a b-d tree */
Int2 BDcompare(PBDNode newnode, PBDNode treeplace, Int2 dimension)
{
	Int2 d;
	FloatLo diff;
PMAD pmad;
PALD pald;
PMGD pmgd;
CharPtr nam;

	for(d=dimension;d<=2;d++) {
		diff=newnode->dim[d]-treeplace->dim[d];
		if (diff!=0.0) return (d+d+(diff>0));
	}
printf(" %f %f %f , %f %f %f\n",newnode->dim[0],newnode->dim[1],newnode->dim[2],treeplace->dim[0],treeplace->dim[1],treeplace->dim[2]);
pald=newnode->paldThis;
pmad=(PMAD)(pald->pfbParent);
pmgd=(PMGD)(pmad->pfbParent);
ErrPostEx(SEV_ERROR,1,1,"..%d..\n",(pmgd->pdnmgLink)->choice);
nam=pmad->pcAName;
ErrPostEx(SEV_ERROR,1,1,"..%s..\n",nam);
ErrPostEx(SEV_ERROR,1,1,"..%d..\n",(((PMGD)(((PMAD)((treeplace->paldThis)->pfbParent))->pfbParent))->pdnmgLink)->choice);
ErrPostEx(SEV_ERROR,1,1,"..%s..\n",((PMAD)((treeplace->paldThis)->pfbParent))->pcAName);
	return(EQUAL);
}

/* inserts newnode into the tree with root treeplace; dimension should be
   zero when called */
void BDInsert(PBDNode newnode,PBDNode treeplace,Int2 dimension)
{
	Int2 branch;
	if (newnode==NULL) return;
	branch=BDcompare(newnode,treeplace,dimension);
	if (branch==EQUAL) {
		ErrPostEx(SEV_FATAL,4,1,"b-d node crash; node not inserted");
		return;
	}
	dimension=branch>>1;
	if(treeplace->way[branch]==NULL) {
		treeplace->way[branch]=newnode;
	}
	else {
		treeplace=treeplace->way[branch];
		BDInsert(newnode,treeplace,dimension);
	}
}

/* this function is mostly for debugging purposes */
void TraverseBDTree(PBDNode pbdTreeHead, Int4 level) {
	static Int4 numnode;
	Int2 cnt;

	if (level==0)
		numnode=0;
	if (pbdTreeHead==NULL)
		return;
	/* non-empty node */
	numnode++;
/*	if (fabs(pbdTreeHead->dim[0]-74.05)<0.005) {
		PrintVec(pbdTreeHead->dim);
		printf("level=%d\n",level);
	}*/
	for (cnt=0;cnt<6;cnt++)
		TraverseBDTree(pbdTreeHead->way[cnt],level+1);
	if (level==0)
		printf("total # nodes: %d\n",numnode);
}

/* remove atoms from Hbond list as they are BDRemoved */
/* if donor or acceptor is BDremoved, no longer a collision so remove from list! */
void ClearHBondStatus(PMAD pmad,PHBS *pphbsHBondsToCheck)
{
	PHBS phbsHere,phbsNext;
	Char atomname;
	
	atomname=(pmad->pcAName)[1];
	/* only affects N and O atoms */
	if (atomname!='N' && atomname!='O')
		return;
	phbsHere=*pphbsHBondsToCheck;
	while (phbsHere!=NULL) {
		phbsNext=phbsHere->next;
		if (pmad==phbsHere->pmadDonor)
			/* clear h-bonds pending */
			(pmad->bReserved)&=0xf1;
		if (pmad==phbsHere->pmadDonor || pmad==phbsHere->pmadAcceptor) {
			if (phbsHere->prev==NULL) {
				*pphbsHBondsToCheck=phbsNext;
				if (*pphbsHBondsToCheck!=NULL)
					(*pphbsHBondsToCheck)->prev=NULL;
			}
			else if (phbsHere->next==NULL) {
				phbsHere->prev->next=phbsNext;
			}		
			else {
				phbsHere->prev->next=phbsNext;
				phbsNext->prev=phbsHere->prev;
			}
			MemFree(phbsHere);
		}
		phbsHere=phbsNext;
	}
}

/* removes killnode->way[which] from the tree with root treeplace */
/* will not remove the root node!!! */
void BDRemove(vec vCoOrd,PBDNode *ppbdTreeHead,PHBS *pphbsHbond)
{
	PBDNode child[6],pbdHere,pbdLast;
	Int2 dimension,cnt,which=0,newhead;
	PALD paldThis;
	PMAD pmadThis;

	pbdHere=*ppbdTreeHead;
	pbdLast=NULL;
	do {
		dimension=0;
		while ((dimension<3) && (pbdHere->dim[dimension]==vCoOrd[dimension]))
			dimension++;
		if (dimension==3) break;
		/* follow right branch */
		pbdLast=pbdHere;
		which=dimension<<1;
		if (vCoOrd[dimension]>pbdHere->dim[dimension]) {
			pbdHere=pbdHere->way[(dimension<<1)+1];
			which++;
		}
		else 
			pbdHere=pbdHere->way[dimension<<1];
	} while (pbdHere);
	if (pbdLast==NULL) {
		/* delete root node */
		if (pbdHere==NULL) return;
		for (cnt=0;cnt<6;cnt++)
			child[cnt]=pbdHere->way[cnt];
		/* clear reserved bit to invalidate co-ordinates at PMAD level */
		paldThis=pbdHere->paldThis;
		pmadThis=(PMAD)(paldThis->pfbParent);
		(pmadThis->bReserved)&=0xfe;
		if (pphbsHbond!=NULL)
			ClearHBondStatus(pmadThis,pphbsHbond);
		MemFree(pbdHere);
		for (cnt=0;cnt<6;cnt++)
			if (child[cnt]!=NULL) 
				break;
		newhead=cnt;
		/* no nodes left? */
		if (newhead==6) {
			*ppbdTreeHead=NULL;		
			return;
		}
		/* re-insert nodes */
		cnt=0;
		/* make way[newhead] the new head */
		pbdLast=child[newhead];
		*ppbdTreeHead=pbdLast;
		while ((cnt<6) && (cnt!=newhead) && (child[cnt]==NULL)) cnt++;
		/* find first non-NULL child if any */
		if (cnt==6)
			return;
		while (cnt<6) {
			/* re-insert remaining children */
			if (child[cnt]!=NULL && cnt!=newhead) {
				BDInsert(child[cnt],pbdLast,0);
			}
			cnt++;
		}
		return;
	}
	if (pbdHere==NULL) {
		ErrPostEx(SEV_ERROR,1,1,"Tried to BDRemove non-existent (%1.5f,%1.5f,%1.5f)",vCoOrd[0],vCoOrd[1],vCoOrd[2]);
		return;
	}
	for (cnt=0;cnt<6;cnt++)
		child[cnt]=pbdHere->way[cnt];
	/* clear reserved bit to invalidate co-ordinates at PMAD level */
	paldThis=pbdHere->paldThis;
	pmadThis=(PMAD)(paldThis->pfbParent);
	(pmadThis->bReserved)&=0xfe;
	if (pphbsHbond!=NULL)
		ClearHBondStatus(pmadThis,pphbsHbond);
	MemFree(pbdHere);
	cnt=0;
	while ((cnt<6) && (child[cnt]==NULL)) cnt++;
	/* find first non-NULL child if any */
	if (cnt==6) {
		pbdLast->way[which]=NULL;
		return;
	}
	pbdLast->way[which]=child[cnt];
	cnt++;
	while (cnt<6) {
		/* re-insert remaining children */
		if (child[cnt]!=NULL) {
			BDInsert(child[cnt],pbdLast->way[which],0);
		}
		cnt++;
	}
}

/* add all atoms in a given PFB and Model to the growing BD-tree with
   pbdTreeHead as its head; if the latter is null, a new tree is created
   and the root of the tree is always returned 
   Adds all conformations (occupanmcy redundancy) if AllConfs is true */
void AddToBDTreeEx(PFB pfbThis, Int2 Model,PBDNode *ppbdTreeHead,Boolean AllConfs)
{
	PMAD pmadAtom;
	PALD paldLoc;
	PBDNode newnode;
	ValNodePtr vnpHead,vnpHere,vnpTail;
	Int4 numnodes;

	if (pfbThis==NULL) return;
	/* get list of atoms */
	vnpHead=ConvertNode(pfbThis,AM_MAD);
	vnpHere=vnpHead;
	if (vnpHere==NULL)
		return;
	/* try adding from the middle of protein and wrap around */
	if (pfbThis->bMe!=AM_MAD) {
		numnodes=1;
		while (vnpHere->next!=NULL) {
			vnpHere=vnpHere->next;
			numnodes++;
		}
		vnpTail=vnpHere;
		vnpHere=vnpHead;
		if (numnodes>1) {
			numnodes=numnodes/2;
			while (numnodes>1) {
				vnpHere=vnpHere->next;
				numnodes--;
			}
			/* now at middle of list, do switch */
			vnpTail->next=vnpHead;
			vnpHead=vnpHere->next;
			vnpHere->next=NULL;
			vnpHere=vnpHead;
		}
	}
	while(vnpHere) {
		/* for each atom in list, get its co-ordinates */
		pmadAtom=(PMAD)(vnpHere->data.ptrvalue);
		paldLoc=GetAtomLocs(pmadAtom,Model);
		while (paldLoc) {
			/* create new node add to tree */
			newnode=(PBDNode)MemNew(sizeof(BDNode));
			newnode->dim[0]=AtomLocX(paldLoc);
			newnode->dim[1]=AtomLocY(paldLoc);
			newnode->dim[2]=AtomLocZ(paldLoc);
			newnode->paldThis=paldLoc;
			/* keep pointer to b-d node in pGraphic */
			paldLoc->pGraphic=newnode;
			/* if root node, set it as root */
			if (*ppbdTreeHead==NULL) {
				*ppbdTreeHead=newnode;
			}
			/* otherwise insert it where it goes */
			else {
				BDInsert(newnode,*ppbdTreeHead,0);
			}
			paldLoc=paldLoc->next;
			if (AllConfs==FALSE)
				break;
		}
		vnpHere=vnpHere->next;
	}
	/* free the temporary atom list we created */
	vnpHead=ValNodeFree(vnpHead);
}

void AddToBDTree(PFB pfbThis, Int2 Model,PBDNode *ppbdTreeHead)
{
	AddToBDTreeEx(pfbThis,Model,ppbdTreeHead,FALSE);
}

/* recursively frees the memory taken by a b-d tree */
void FreeBDTree(PBDNode pbdTreeHead)
{
	PBDNode next[6];
	Int2 cnt;

	if (pbdTreeHead==NULL) return;
	for (cnt=0;cnt<6;cnt++)
		next[cnt]=pbdTreeHead->way[cnt];
	pbdTreeHead=MemFree(pbdTreeHead);
	/* recursively free the b-d Tree */
	for (cnt=0;cnt<6;cnt++)
		FreeBDTree(next[cnt]);				
}

/* adds a new world pointed to by pwsNew to the current world list, headed
   by global variable WorldsHead, returning a pointer to it; it is
   assigned an Id fo later retrieval with FindWorld */
ValNodePtr InstantiateWorld(Int2 Id,PWS pwsNew)
{
	ValNodePtr NewWorld;

	/* add the world to the linked list of worlds */
	NewWorld=ValNodeAddPointer(&WorldsHead,Id,pwsNew);
	return NewWorld;
}

/* frees the memory used by a given world and its associated
   sub-structures; the argument should be the value returned by
   InstantiateWorld given when that world was instantiated; the return
   value is always NULL */
ValNodePtr FreeWorld(ValNodePtr vnpKill)
{
        ValNodePtr vnpHere;
        PWS pwsHead,pwsNext;        
   
	if (vnpKill==NULL) return NULL;
	vnpHere=WorldsHead;
	/* unlink the world node */
	if (vnpHere==vnpKill)
		WorldsHead=vnpHere->next;
	else {
		while(vnpHere->next != vnpKill)
			vnpHere=vnpHere->next;
		vnpHere->next=vnpKill->next;
	}
	pwsHead=(PWS)(vnpKill->data.ptrvalue);
	/* free the single Val Node */
	vnpKill=MemFree(vnpKill);
	if (pwsHead==NULL) return NULL;
	/* free b-d tree hanging off first node */
	FreeBDTree(pwsHead->pbdTreeHead);
        /* walk down list deleting one node at a time, including molecule
	   list if it exists */
        while (pwsHead) {
                pwsNext=pwsHead->next;
                MemFree(pwsHead);
                pwsHead=pwsNext;
        }
        return NULL;   
}

/* frees the memory used by all worlds and associated sub-structures; it
   should be called at the end of any program which has used worlds */
ValNodePtr FreeAllWorlds(void)
{
	while(WorldsHead)
		FreeWorld(WorldsHead);
	return NULL;
}

/* if pwsCurrent is NULL, a new linked list is formed, otherwise the node
   is appended to the end of pwsCurrent; the new node is
   always returned.  Adds a PFB in a given Model to the list of objects in
   this world; the PFB must be a PMSD, PMMD, PMGD, PMAD or PMBD or the
   results will be unpredictable (and incorrect) */
PWS AddtoWorldEx(PWS pwsCurrent,Int2 Model,PFB pfbThis,Boolean AllConfs)
{
	PWS pwsNode,pwsHere;

	/* create the new world node and set its fields */
	pwsNode=(PWS)MemNew(sizeof(WS));
	pwsNode->next=NULL;
	/* link new node to end of current node list */
	if (pwsCurrent!=NULL) {
		pwsHere=pwsCurrent;
		while(pwsHere->next)
			pwsHere=pwsHere->next;
		pwsHere->next=pwsNode;
	}
	pwsNode->Model=Model;
	pwsNode->pfbThis=pfbThis;
	pwsNode->pbdTreeHead=NULL;
	/* add new points in this node to b-d tree, which is hanging off
	   this node */
	/* return the head of the list, which is unchanged unless we just
	   created the first node */
	if (pwsCurrent==NULL)
		AddToBDTreeEx(pwsNode->pfbThis,Model,&(pwsNode->pbdTreeHead),AllConfs);
	else
		AddToBDTreeEx(pwsNode->pfbThis,Model,&(pwsCurrent->pbdTreeHead),AllConfs);
	return pwsNode;
}

PWS AddtoWorld(PWS pwsCurrent,Int2 Model,PFB pfbThis)
{
	return AddtoWorldEx(pwsCurrent,Model,pfbThis,FALSE);
}
/* locates and returns a pointer to the first instantiated world
   with given ID */
PWS FindWorld(Int2 ID)
{
	ValNodePtr vnpHead;

	vnpHead=WorldsHead;
	/* walk through list comparing IDs until we find the right one */
	while (vnpHead) {
		if (vnpHead->choice==ID)
			return (PWS)(vnpHead->data.ptrvalue); 
		vnpHead=vnpHead->next;
	}
	return NULL;
}

/* used recursively, this traverses the tree with root pbdTreeplace and
   returns a pointer to a list of ValNodes, pointed to by vnptarget,
   containing PALDs whose centres are within a sphere centred at
   (min+max)/2 and with radius (max[0]-min[0])/2; dimension should be 0
   when called; if vnptarget is NULL, the head of the new list is
   returned, otherwise new nodes are appended to the list pointed to by
   vnptarget, and the return value is the same vnptarget */ 
void DimNeighbor(PBDNode pbdTreeplace,Int2 dimension)
{
	static PBDNode pbdTemp;
	static FloatLo distsq,dist1;
	static FloatLo *coord;
	
	if (pbdTreeplace==NULL) return;
	/* check if left branch valid */
	coord=pbdTreeplace->dim;
	if (coord[dimension]>minbound[dimension])
		if ((pbdTemp=pbdTreeplace->way[dimension<<1])!=NULL)
			DimNeighbor(pbdTemp,dimension);
	/* recurse to one dimension deeper */
	if (dimension<2)
		DimNeighbor(pbdTreeplace,dimension+1);
	/* if reach here, could be a valid point at this node */
	else if ((coord[0]<maxbound[0]) && (coord[1]<maxbound[1]) && (coord[2]<maxbound[2])) {
		/* calculate distance and compare to given parameters */
		dist1=coord[0]-sumbound[0];
		distsq=dist1*dist1;
		dist1=coord[1]-sumbound[1];
		distsq+=(dist1*dist1);
		dist1=coord[2]-sumbound[2];
		distsq+=(dist1*dist1);
		/* and add to list if inside the bounding sphere */
		if (distsq<circradsq) {
			vnpNearAtomsLast=ValNodeNew(vnpNearAtomsLast);
			vnpNearAtomsLast->data.ptrvalue=(void *)(pbdTreeplace->paldThis);
/*			ValNodeAddPointer(&vnpNearAtoms,0,pbdTreeplace->paldThis);
*/		}
	}
	/* check if right branch valid */
	if (coord[dimension]<maxbound[dimension])
		if ((pbdTemp=pbdTreeplace->way[1+(dimension<<1)])!=NULL)
			DimNeighbor(pbdTemp,dimension);
	return;
}

/* finds all atoms in the world pwsWorld whose centres lie in a sphere of
   given radius and centre, returning them as a ValNodePtr to a list of
   PALDs */
ValNodePtr FindAtomsIn(PWS pwsWorld,vec centre,FloatLo radius)
{
	PBDNode pbdTreeHead;
	Int2 cnt;

	if (pwsWorld==NULL) return NULL;
	pbdTreeHead=pwsWorld->pbdTreeHead;
	if (pbdTreeHead==NULL) return NULL;
	/* find limits of a bounding box needed by DimNeighbor, store in static globals */
	for (cnt=0;cnt<3;cnt++) {
		minbound[cnt]=centre[cnt]-radius;
		maxbound[cnt]=centre[cnt]+radius;
	}
	/* traverse b-d tree adding ValNodes for each valid PALD */
	/* sentinel node at head of list */
	vnpNearAtoms=ValNodeNew(NULL);
	vnpNearAtomsLast=vnpNearAtoms;
	/* precalculate constants */	
	circradsq=(maxbound[0]-minbound[0])/2.0;
	circradsq=circradsq*circradsq;
	sumbound[0]=(maxbound[0]+minbound[0])/2.0;
	sumbound[1]=(maxbound[1]+minbound[1])/2.0;
	sumbound[2]=(maxbound[2]+minbound[2])/2.0;
	DimNeighbor(pbdTreeHead,0);
	/* free sentinel node */
	vnpNearAtomsLast=vnpNearAtoms;
	vnpNearAtoms=vnpNearAtoms->next;
	MemFree(vnpNearAtomsLast);
	return vnpNearAtoms;
}
/*
void DumNeighbor(PBDNode pbdTreeplace,Int2 dimension)
{
	Int2 cnt;
	FloatLo distsq,dist1;

	if (pbdTreeplace==NULL) return;
	pbdHere=pbdTreeHead;
	while (pbdHere->dim[dimension]>minbound[dimension] && pbdHere->way[dimension<<1]!=NULL)
		pbdHere=pbdHere->way[dimension<<1];
	if (dimension<2)
		DumNeighbor(pbdHere,dimension+1);
	else if ((coord[0]<maxbound[0]) && (coord[1]<maxbound[1]) && (coord[2]<maxbound[2])) {
		distsq=0;
		for (cnt=0;cnt<3;cnt++) {
			dist1=coord[cnt]-sumbound[cnt];
			distsq=distsq+dist1*dist1;
		}
		if (distsq<circradsq) {
			ValNodeAddPointer(&vnpNearAtoms,0,pbdTreeplace->paldThis);
		}
	}
	if (coord[dimension]<maxbound[dimension])
		DimNeighbor(pbdTreeplace->way[1+(dimension<<1)],dimension);
	return;
}
*/


/*  
$Log: bbox.c,v $
Revision 1.20  2003/12/08 22:52:27  hfeldman
Fixed typo

Revision 1.19  2003/12/08 19:09:25  hfeldman
Fixed potential crash

Revision 1.18  2003/11/06 22:29:33  feldman
Fixed compiler warning

Revision 1.17  2003/10/01 17:04:38  feldman
Added occupancy redundancy as optional flag for adding to BD tree

Revision 1.16  2001/11/02 20:32:44  feldman
Fixed compiler warnings

Revision 1.15  2001/10/23 15:25:41  feldman
Fixed h-bonding bug

Revision 1.14  2001/10/10 19:48:19  feldman
Revamped hydrogen bonding code - now much more robust and most likely correct

Revision 1.13  2001/10/02 21:52:51  feldman
commented out error msg

Revision 1.12  2001/10/01 14:46:02  feldman
Changed PRECISION to 0.0 for BD tree

Revision 1.11  2001/09/26 15:31:26  feldman
Added Sovlation term to energy for EEF1

Revision 1.10  2001/07/26 19:19:51  feldman
Changed way fragments are stored on disk to be more platform independent,
and removed id parameter from AddtoWorld

Revision 1.9  2001/07/13 19:01:14  feldman
Id field removed from bd functions

Revision 1.8  2001/07/12 21:06:48  feldman
started changing to allow deletion of root node

Revision 1.7  2001/03/14 16:25:55  feldman
Updated headers, adding license where missing, changing &quot to real
quote, added extern C where needed and #ifdef FILENAME_H etc.,
plus got rid of some unneeded variables/functions

Revision 1.6  2000/09/26 18:33:07  feldman
removed endless loop

Revision 1.5  2000/07/13 16:05:16  feldman
Removed further dependencies between foldtraj and general structure
functions

Removed all FATAL errors from non-foldtraj c files (except one)
and PurgeGlobs
Most programs should compile without needing libfoldtraj.a now

Revision 1.4  2000/07/07 21:30:27  feldman
Tidied up .h header files and removed some
inter-dependencies

Revision 1.3  2000/06/20 16:40:23  feldman
Incorporated sstru extended structure calculation into C code,
no longer need an external program

Revision 1.2  2000/06/15 15:52:10  feldman
Corrected Bzip2 log dumping and improved H-bond treatment (fixed
bug that was crashing program)

Revision 1.1.1.1  2000/06/09 18:13:58  feldman
TraDES project

*/

