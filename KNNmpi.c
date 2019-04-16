/**************************************************************************************************************************************************************************************
/************************************************************************************************************************************ ***** */
/********************************************************************************************************************************************************************************* */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#include <sys/time.h>
#include <stddef.h>
#define SIZE0 3 // stoixeia dianismatos me arg[v]. meta.....
#define INIT_MINDIST 1000000

typedef struct Points {
    double x;
    float y;
	float z;
	float maxDist;
	int type; //Nq or Nc
	int bx;
	int by;
	int bz;
	int process;
	float nbx;float nby;float nbz;
}point;

int n,m,k;
int rc,numTasks,rank;
time_t timee	;

//FUNCTIONS 
void sendToproc(int *forSend,int pntx,int pnty,int pntz,int nmb,int col,int row,int *sendCout);
int getProc(int bx,int by,int bz,int pnum,int n,int m,int k);
float getDistance(float x2,float y2,float z2,point p);
void measure(point *arr,point *a[n][m][k],int counter[n][m][k],int n,int m,int k,int j,int row,int column,int *sendCout,int numtasks);
void measure2(point *arr,point *a[n][m][k],int counter[n][m][k],int px,int py,int pz,int j,int row,int column,int *sendCout,int numtasks,int n,int m,int k);
float *neighbours(point *arr,int counter[n][m][k],int n,int m,int k,int pnum,int *sendCout,int column,float *forSend);
/************************************************************ STARTING MAIN PROGRAM **************************************************************************************************/
int main(int argc,char *argv[]){
	
int Nc,Nq;
struct timeval start0,end0,start1,end1,start2,end2;
int j,n,m,k,numTasks,SelfTID;
if (argc != 6) {
printf("Usage: %s q\n  where n=2^q is problem size (power of two)\n", 
argv[0]);
exit(1);
}	
n = 1<<atoi(argv[1]);
m = 1<<atoi(argv[2]);
k = 1<<atoi(argv[3]); 

Nc=1<<atoi(argv[4]); 
Nq=1<<atoi(argv[5]); 

MPI_Status status[4];
MPI_Request req[4];
rc = MPI_Init(&argc,&argv);
if (rc != MPI_SUCCESS) {
printf ("Error starting MPI program. Terminating.\n");
MPI_Abort(MPI_COMM_WORLD, rc);
}
MPI_Comm_size(MPI_COMM_WORLD,&numTasks);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);

srand((unsigned) time(&timee)+200*rank);  
float tmp;
int i,cc,cout,rc;
float tx,ty,tz;  //Grid's boxes coordinates
int size=((Nc+Nq)/2)-1;
int	pointNum=(Nc+Nq)/numTasks;

point *arr=(point *)malloc((Nc+Nq+128)*sizeof(point));   //Allocate memory for local array holding points from current process

//set number 1	1st half=Q type of points,2nd half = C type of points
point *a[n][m][k];
int counter[n][m][k];
for(i=0;i<n;i++){
  for(j=0;j<m;j++){
   for(cc=0;cc<k;cc++){
	a[i][j][cc]=(point *)malloc((Nc+Nq/n*m*k +15)*sizeof(point));
	counter[i][j][cc]=0;
}
}
}
	
int sum=0;

for(i=0;i<pointNum;i++){

  sum+=1;
  tmp=((float)rand())/(RAND_MAX);
  arr[i].x=tmp;
  tmp=((float)rand())/(RAND_MAX);        //Store random values (uniform distribution) for 3 coordinates of each point
  arr[i].y=tmp;
  tmp=((float)rand())/(RAND_MAX);
  arr[i].z=tmp;

  if(i<=(pointNum/2))
  	arr[i].type=1;
  else
  	arr[i].type=2;
 
  float n1=(float)n;
  float m1=(float)m;
  float k1=(float)k;

  
  //Store each point in box that belongs (in unit cube)
  //x coordinate
if(fmod((arr[i].x),(1.0/n1))==0){
	if(arr[i].x<(1.0/n1)){
	tx=0;
}
else{
	tx=(arr[i].x/(1.0/n1))-1;	
}
}
else{
	tx=(arr[i].x/(1.0/n1));
}
//y coordinate
if((fmod(arr[i].y,(1.0/m1)))==0){
if((arr[i].y)<(1.0/m1)){
	ty=0;
}
else{
	ty=((arr[i].y)/(1.0/m1))-1;
}	
//(arr[i].y)fmod(pow(2,-m)==0			
}
else{
	ty=((arr[i].y)/(1.0/m1));
}
//z coordinate
if(fmod(arr[i].z,(1.0/k1))==0){
if((arr[i].z)<(1.0/k1)){
	tz=0;
}
else{
	tz=((arr[i].z)/(1.0/k1))-1;			
}	
}
else{
	tz=((arr[i].z)/(1.0/k1));
} 
int t1=(int)tx; //coordinates of box(in grid) that current
int t2=(int)ty; //point belongs
int t3=(int)tz;  //increasing counter(frequency..)

//Passin value of points in box' array of structs 
arr[i].bx=t1;
arr[i].by=t2;
arr[i].bz=t3;
arr[i].process=getProc(t1,t2,t3,numTasks,n,m,k);
arr[i].maxDist=INIT_MINDIST;
cout=counter[t1][t2][t3];

((a[t1][t2][t3])[cout])=arr[i];
((a[t1][t2][t3])[cout]).x=arr[i].x;
((a[t1][t2][t3])[cout]).y=arr[i].y;
((a[t1][t2][t3])[cout]).z=arr[i].z;
((a[t1][t2][t3])[cout]).bx=arr[i].bx;
((a[t1][t2][t3])[cout]).by=arr[i].by;
((a[t1][t2][t3])[cout]).bz=arr[i].bz;
((a[t1][t2][t3])[cout]).process=arr[i].process;
counter[t1][t2][t3]+=1;

((a[t1][t2][t3])[cout]).process=arr[i].process;
}


int max1=-2;
for(i=0;i<n;i++){
 for(j=0;j<m;j++){
 for(cc=0;cc<k;cc++){
	if(counter[i][j][cc]>max1)  
		max1=counter[i][j][cc];
}}}
/*create struct type for MPI */
const int items=9;
int 		   blocklengths[9]={1,1,1,1,1,1,1,1,1};
MPI_Datatype   types[12]={MPI_FLOAT,MPI_FLOAT,MPI_FLOAT,MPI_FLOAT,MPI_INT,MPI_INT,MPI_INT,MPI_INT,MPI_INT,MPI_FLOAT,MPI_FLOAT,MPI_FLOAT};
MPI_Datatype   mpi_point_type;
MPI_Aint       offsets[12];

offsets[0]=offsetof(point,x);
offsets[1]=offsetof(point,y);
offsets[2]=offsetof(point,z);
offsets[3]=offsetof(point,maxDist);
offsets[4]=offsetof(point,type);
offsets[5]=offsetof(point,bx);
offsets[6]=offsetof(point,by);
offsets[7]=offsetof(point,bz);
offsets[8]=offsetof(point,process);
offsets[9]=offsetof(point,nbx);
offsets[10]=offsetof(point,nby);
offsets[11]=offsetof(point,nbz);
MPI_Type_create_struct(items,blocklengths,offsets,types,&mpi_point_type);
MPI_Type_commit(&mpi_point_type);
  //MPI_Type_indexed(3, blocklengths, o, type2, &type1);
  //  MPI_Type_commit(&type1);
point recBuf,tBuf;
int cout0=0,tag=221,tmpCout;  //tmpCout holds new number of cells for each grid's box
int cox,coy,coz,flag,newCout;

//****************************************************************************************************************************
//****************************************************************************************************************************
//****************************************************************************************************************************
//Sending and receiving points from one process to another

for(i=0;i<pointNum;i++){
  int proc;
  proc=getProc(arr[i].bx,arr[i].by,arr[i].bz,numTasks,n,m,k);//Get Process that each point(of current process) belongs..
  arr[i].process=proc;  //get process number for each point of current process
  
  if(proc!=rank){    //If point of this process dont match
  	point tm=arr[i];	//area criteria for current process..Send to right process
  	MPI_Isend(&tm,1,mpi_point_type,proc,0,MPI_COMM_WORLD,&req[0]);
  }
  
  while(1){
  	MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&flag,&status[0]);
  	if(!flag) break; //if message had not arrived..break from  
  				//infin. loop
  	MPI_Irecv(&recBuf,1,mpi_point_type,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&req[0]);
  	arr[pointNum+cout0]=recBuf;
   
     //ERROR CHECKING
    int errChk =getProc(recBuf.bx,recBuf.by,recBuf.bz,numTasks,n,m,k);
    if(errChk!=rank){
      printf("ERROR RECEIVED POINT DOESNT BELONG TO CURRENT PROCESS \n");

    }
    
  	cout0+=1; 					//increase counter of array after each point receive
  	tmpCout=counter[recBuf.bx][recBuf.by][recBuf.bz];
  	((a[recBuf.bx][recBuf.by][recBuf.bz])[tmpCout])=recBuf;
  	counter[recBuf.bx][recBuf.by][recBuf.bz]+=1;
  	//if(rank==0)
  	//	printf("ADDED in box %d \n",(((a[recBuf.bx][recBuf.by][recBuf.bz])[tmpCout]).bz));
  
  }
}//end for
//*******************************************************************************new code*******************************************************************************
//****************************************************************************************************************************
//****************************************************************************************************************************
//****************************************************************************************************************************

MPI_Barrier(MPI_COMM_WORLD);
 
 //Start timer  , in rank0
 if(rank==0)
    gettimeofday(&start0, NULL);
    
if(rank==0) printf("HELLO\n");
point *res;
float fdist;
int row=numTasks;int totalNum;

int col=Nq+256;  //Size for possible candidates to send in each other process

int *sendCout=malloc(numTasks*sizeof(int));
for(i=0;i<numTasks;i++){
	sendCout[i]=0;}
float *forSend;
if((forSend=(float *)malloc(numTasks*col*sizeof(float)))==NULL){ // 2d array stored in 1d (points sent between processes) ,1 neighbour..(3 FLOATS PER POINT)!!!!
	printf("NOT ENOUGH MEMORY \n");
	exit(1);
}

float *rcvbuf;
if((rcvbuf=(float *)malloc(numTasks*col*sizeof(float)))==NULL){ // 2d array stored in 1d (points sent between processes) ,1 neighbour..(3 FLOATS PER POINT)!!!!
	printf("NOT ENOUGH MEMORY \n");
	exit(1);
}


totalNum=pointNum+cout0;

//Function for finding neighbour processes and sending slices of points 
rcvbuf=neighbours(arr,counter,n,m,k,totalNum,sendCout,col,forSend);  
int coutt=0;
for(j=0;j<(totalNum);j++){
	if((arr[j].type=1)&&(arr[j].process==rank)){	
		
	//First checking points that contains our local array
	measure(arr,a,counter,n,m,k,j,row,col,sendCout,numTasks);//Function for finding (LOCAL)  min distance(min distance return)
	
	
	//After checking points received from neighbours ( IN rcvbuf)
	if(numTasks<32){
		for(i=0;i<col;i=i+3){
			if((rank-1)>0){
				fdist=getDistance(rcvbuf[(rank-1)*col+i],forSend[(rank-1)*col+i+1],rcvbuf[(rank-1)*col+(i+2)],arr[j]);
				if((fdist<(arr[j].maxDist))&&(fdist!=0))    // !! Ensure not to check current point with its self (!)
				arr[j].maxDist=fdist;
			}
			if((rank+1)<numTasks){
				fdist=getDistance(rcvbuf[(rank+1)*col+i],rcvbuf[(rank+1)*col+(i+1)],rcvbuf[(rank+1)*col+(i+2)],arr[j]);
				if((fdist<(arr[j].maxDist))&&(fdist!=0))  
				arr[j].maxDist=fdist;
			}
		}
	}
	
 
  //RESULT PRINTING (FOR RANK=0) .MANUAL CHANGE FOR PRINTING FOR ALL PROCESSES
 	if(rank==0){
	if(arr[j].maxDist!=INIT_MINDIST)
	printf("FOR RANK : %d , MIN.DISTANCE OF POINT %f %f %f  is %f  .NEIGHBOUR IS : %f %f %f \n",rank,arr[j].x,arr[j].y,arr[j].z,arr[j].maxDist,arr[j].nbx,arr[j].nby,arr[j].nbz	);
	}
	
	
	}
}
MPI_Barrier(MPI_COMM_WORLD);
if(rank==0)
printf("FINISHED\n");

//***********************start new code send-receive for measure


if(rank==0)  {
gettimeofday(&end0, NULL);
printf("FOR FIRST KERNEL, VERSION WITHOUT SHARED MEMORY TIME IS : %ld microSeconds \n", ((end0.tv_sec * 1000000 + end0.tv_usec)
		  - (start0.tv_sec * 1000000 + start0.tv_usec)));
}

free(arr);
//free(recvBuf);
free(sendCout);
MPI_Type_free(&mpi_point_type);

MPI_Finalize();

return 0;
} //end MAIN*********************************************************************************************************************************************************************
//*******************************************************************************************************************************************************************************
float getDistance(float x2,float y2,float z2,point p){
float distance;

float tx=(p.x)-x2;
float tx0=pow(tx,2);

float ty=(p.y)-y2;
float ty0=pow(ty,2);

float tz=(p.z)-z2;
float tz0=pow(tz,2);

distance=sqrt(tx0+ty0+tz0);
return(distance);

}
//********************************************************************************************************s
//cx,cy,cx ->bx,by,bz (box coordinates)
int getProc(int bx,int by,int bz,int pnum,int n,int m,int k){
int i,j,l;
int temp;
if(pnum==2){
if(bz<(k/2))
return 0;   //process 0 ->n,m whole interval,z->0,k-1
else
return 1;
}

if((pnum>2)&&(pnum<32)){  //minimum value of k(grid's coord) is 2^4 so at least
temp=k/pnum; //Distance between boxes in grid's k dimension (z axis)
//one plane(in z axis)will be available for each process
return(bz/temp);  // Process that given point belongs
}
//	*****************************************************
if(pnum==32){
temp=k/16;   //distance between process area,in z axis
if(by>(m/2))	//in X-axis,process cut the area  (16,16)
return((pnum/2)+((bz)/temp)); //Process number` start counting from zero (0)
else
return((bz)/temp);
}
if(pnum==64){
temp=k/16;
if((by)<(m/4))
	return((bz)/temp);			//Divides grid in 4 equal pieces
else if((by)<(m/2))					//in y axis(m size) of grid
	return(((bz)/temp)+(pnum/4));	//and as before in z axis(k size of grid)
else if((by)<(0.75*m))
	return(((bz)/temp)+(pnum/2));
else
	return(((bz)/temp)+(0.75*pnum));
}										

}
//********************************************************************************************************

//function that calculates neighbours and minimum distance for current point given 
//It looks in same box of grid,  and in ALL 27 neighbour boxes

void measure(point *arr,point *a[n][m][k],int counter[n][m][k],int n,int m,int k,int j,int column,int row,int *sendCout,int numtasks){
int flag2,numb,tmpx,tmpy,tmpz=0,i;
float dist;
int procPos;  //process of box to be sent
//Allocation of 2d array that stores points to be sent to other processes

int cx=arr[j].bx;
int cy=arr[j].by;				
int cz=arr[j].bz;int num;
 //if((rank==0)&&(cz>k/2+2))
 //printf(" getting point %d %d %d from main \n",cx,cy,cz);
  
int maxC=counter[cx][cy][cz]; //New counter after point exchanges between processess 
for(i=0;i<maxC;i++){
	//searching in current box...	
	if(((a[cx][cy][cz])[i]).type==2){
	int zz= ((a[cx][cy][cz])[i]).type;
	dist=getDistance((((a[cx][cy][cz])[i]).x),(((a[cx][cy][cz])[i]).y),(((a[cx][cy][cz])[i]).z),arr[j]);
	if((dist<(arr[j].maxDist))&&(dist!=0)) { //minDist..
		arr[j].nbx=((a[cx][cy][cz])[i]).x;arr[j].nby=((a[cx][cy][cz])[i]).y;arr[j].nbz=((a[cx][cy][cz])[i]).z;
		arr[j].maxDist=dist;
	}
	}
}

	//CHECKING if our box is in-bounds . If yes ,look for neighbours in those neighour boxes
	
 if((cx>=1) || (cx<n-1)){
	if(cx>=1){		
			measure2(arr,a,counter,(cx-1),cy,cz,j,column,row,sendCout,numtasks,n,m,k);
	if(cy>0){
			measure2(arr,a,counter,(cx-1),(cy-1),cz,j,column,row,sendCout,numtasks,n,m,k);
	if(cz>0){
			measure2(arr,a,counter,(cx-1),(cy-1),(cz-1),j,column,row,sendCout,numtasks,n,m,k);
	}//cz>0
	if(cz<k-1){
			measure2(arr,a,counter,(cx-1),(cy-1),(cz+1),j,column,row,sendCout,numtasks,n,m,k);
	}//cz<k-1
	
	}//	cy>0
	
	if(cy<m-1){
			measure2(arr,a,counter,(cx-1),(cy+1),cz,j,column,row,sendCout,numtasks,n,m,k);
	if(cz>0){
			measure2(arr,a,counter,(cx-1),(cy+1),(cz-1),j,column,row,sendCout,numtasks,n,m,k);
	}
	if(cz<k-1){
			measure2(arr,a,counter,(cx-1),(cy+1),(cz+1),j,column,row,sendCout,numtasks,n,m,k);	
	}
	
	}//cy<m-1
	if(cz>0){
			measure2(arr,a,counter,(cx-1),cy,(cz-1),j,column,row,sendCout,numtasks,n,m,k);
	}
	if(cz<k-1){
			measure2(arr,a,counter,(cx-1),cy,(cz+1),j,column,row,sendCout,numtasks,n,m,k);
	}
	}//cx>=1
	if(cx<n-1){
			measure2(arr,a,counter,(cx+1),cy,cz,j,column,row,sendCout,numtasks,n,m,k);
	if(cz<k-1){
			measure2(arr,a,counter,(cx+1),cy,(cz+1),j,column,row,sendCout,numtasks,n,m,k);
			if(cy<m-1)
				measure2(arr,a,counter,(cx+1),(cy+1),(cz+1),j,column,row,sendCout,numtasks,n,m,k); //NEW
	}
	if(cz>0){
			measure2(arr,a,counter,(cx+1),cy,(cz-1),j,column,row,sendCout,numtasks,n,m,k);
	}
	if(cy>0){
			measure2(arr,a,counter,(cx+1),(cy-1),cz,j,column,row,sendCout,numtasks,n,m,k);	
			if(cz>0){
			measure2(arr,a,counter,(cx+1),(cy-1),(cz-1),j,column,row,sendCout,numtasks,n,m,k);	
			}
		if(cz<k-1)
			measure2(arr,a,counter,(cx+1),(cy-1),(cz+1),j,column,row,sendCout,numtasks,n,m,k);	
	}
	if(cy<m-1){
			measure2(arr,a,counter,cx+1,cy+1,cz,j,column,row,sendCout,numtasks,n,m,k);	
	}
	
	
	}//cx<n-1
	if(cy<m-1){
			measure2(arr,a,counter,cx,cy+1,cz,j,column,row,sendCout,numtasks,n,m,k);
			if(cz>0)
				measure2(arr,a,counter,cx,cy+1,cz-1,j,column,row,sendCout,numtasks,n,m,k);	
			if(cz<k-1)
				measure2(arr,a,counter,cx,cy+1,cz+1,j,column,row,sendCout,numtasks,n,m,k);	
	}
	if(cy>0){
			measure2(arr,a,counter,cx,cy-1,cz,j,column,row,sendCout,numtasks,n,m,k);
			if(cz>0)
			measure2(arr,a,counter,cx,cy-1,cz-1,j,column,row,sendCout,numtasks,n,m,k);
			if(cz<k-1)
			measure2(arr,a,counter,cx,cy-1,cz+1,j,column,row,sendCout,numtasks,n,m,k);
	}
	if(cz>0){
			measure2(arr,a,counter,cx,cy-1,cz,j,column,row,sendCout,numtasks,n,m,k);
			measure2(arr,a,counter,cx,cy+1,cz,j,column,row,sendCout,numtasks,n,m,k);
			}
	
 }//cx>=1||<n-1



}//end measure

//This function computes distance between our current point running in loop (main code) and those from neighbours' boxes 
// !!!					Called from measure function

void measure2(point *arr,point *a[n][m][k],int counter[n][m][k],int px,int py,int pz,int j,int column,int row,int *sendCout,int numtasks,int n,int m,int k){
	int numb,i,procPos;
	float dist;
	if(getProc(px,py,pz,numtasks,n,m,k)==rank){    //if point belongs in other process : Do nothing
			numb=counter[px][py][pz];
			for(i=0;i<numb;i++){
			if((((a[px][py][pz])[i]).type==2)){
			dist=getDistance(((a[px][py][pz])[i]).x,((a[px][py][pz])[i]).y	,(((a[px][py][pz])[i]).z),arr[j]);
			//	printf("2nd Distance %f..\n",dist);
			if((dist<(arr[j].maxDist))&&(dist!=0)) {    // !! Ensure that its not checking point with its self (!)
				arr[j].maxDist=dist;
				arr[j].nbx=((a[px][py][pz])[i]).x;arr[j].nby=((a[px][py][pz])[i]).y;arr[j].nbz=((a[px][py][pz])[i]).z;
			}
			}
		}
		//If points belongs to current process
		//******************************************************************
		
		
	//	
	//pp=getProc(px,py,pz.proc) . Go to array[pp][i]...Dist
			
	}
	
}

//This function finds neighbour processes and send them candidates for nearest neighbours (slices of points that have contact with those processes
//Also receives as described above

float *neighbours(point *arr,int counter[n][m][k],int n,int m,int k,int pnum,int *sendCout,int column,float *forSend){	
	MPI_Request request,request0;
	MPI_Status stats;
	
	int numtasks,i,j,jj,pntx,pnty,pntz,cout0,z0,pos,num;//Cout0 -Counter of current points in forSend array
	float coordinates;
	MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
 	float *recvBuf=(float*)malloc(column*numtasks*sizeof(float));  

	for(i=0;i<pnum;i++){
		if(arr[i].type=2){
			pntx=arr[i].bx;pnty=arr[i].by;pntz=arr[i].bz;	
			if(numtasks==2){

					if((rank==1)&&(pntz==k/2)&&(getProc(pntx,pnty,pntz-1,numtasks,n,m,k)!=rank)){   //eimai sto 2o process se synoriako shmeio..Stelnw (rank-1)
					num=sendCout[rank-1];
					   //Store coordinates of neighbour point
					forSend[((rank-1)*column)+num]=arr[i].x;forSend[(rank-1)*column+(num+1)]=arr[i].y;forSend[(rank-1)*column+(num+2)]=arr[i].z;			//in 32 bit register (10 bit for each coordinate in axis(x,y,z)
			//		printf("%f in box %d ...\n",forSend[(rank-1)*column+(num+2)],pntz);
					sendCout[0]+=3;
					}
					if((rank==0)&&(pntz==k/2-1)&&(getProc(pntx,pnty,pntz+1,numtasks,n,m,k)!=rank)){   //eimai sto 1o process.Se synoriako shmeio stelnw sto rank+1
					num=sendCout[1];
					forSend[column+num]=arr[i].x;forSend[(column+(num+1))]=arr[i].y;forSend[column+(num+2)]=arr[i].z;
				//	printf("%f in box %d ... IN PLACE %d .COLUMN=%d\n",forSend[column+(num+2)],pntz,column+num,column);
					sendCout[1]+=3;
					}	
			} //endif numtasks=2
		
		
		}
	}
		
	
	
		MPI_Datatype rowtype;
		MPI_Type_contiguous(column,MPI_FLOAT,&rowtype);
		MPI_Type_commit(&rowtype);		
		MPI_Barrier(MPI_COMM_WORLD);
		for(j=0;j<numtasks;j++){
			if(j!=rank){
				MPI_Isend(&forSend[j*column],1,rowtype,j,0,MPI_COMM_WORLD,&request);
				MPI_Irecv((recvBuf+column*j),1,rowtype,j,MPI_ANY_TAG,MPI_COMM_WORLD,&request);	//Receive starts ,where sending process' address 
				//test printf (ensure that points have been received)
		//		if(rank==0){	
		//		for(i=0;i<250;i++){	
		//			printf("%f bla bla\n",recvBuf[1*col+i]);				
		//		} 
				
			}//end if j!=rank

		}//end for ...numTasks	
   		//		MPI_Barrier(MPI_COMM_WORLD);
//		MPI_Type_free(&rowtype);
	return recvBuf;
}
//**********************************************************************************************************
void sendToproc(int *forSend,int pntx,int pnty,int pntz,int position,int col,int row,int *sendCout){

	int num=sendCout[position];
	if(rank==0)
		printf("NUM=%d ---rank= %d--- position=%d  point : %d , %d , %d\n",sendCout[position],rank,position,pntx,pnty,pntz);
		
	float coordinates = pntx | (pnty << 10) | (pntz << 20);
	forSend[position*col+num]=coordinates;
	
	sendCout[position]+=1;
}
//**********************************************************************************************************
