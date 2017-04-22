#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stddef.h>
#define SIZE0 3 // stoixeia dianismatos me arg[v]. meta.....
#define INIT_MINDIST 1000000
//trial version of current point recognition 
typedef struct Area{
	float x1[2];
	float y1[2];
	float z1[2];
}areas;
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

typedef struct Boxes{
	float coord[3];
	int cell;	
}box;

int n,m,k;
int SIZE1,SIZE2;
int rc,numTasks,rank;
time_t timee	;
//NEWWW
void sendToproc(int *forSend,int pntx,int pnty,int pntz,int nmb,int col,int row,int *sendCout);
int getProc(int bx,int by,int bz,int pnum,int n,int m,int k);
float getDistance(float x2,float y2,float z2,point p);
void measure(point *arr,point *a[n][m][k],int counter[n][m][k],int n,int m,int k,int j,int row,int column,int *sendCout,int numtasks);
void measure2(point *arr,point *a[n][m][k],int counter[n][m][k],int px,int py,int pz,int j,int row,int column,int *sendCout,int numtasks,int n,int m,int k);

int main(int argc,char *argv[]){
	
struct timeval start0,end0,start1,end1,start2,end2;
int j,n,m,k,numTasks,SelfTID;
if (argc !=6) {
printf("Usage: %s \n  where n=2^q is problem size (power of two)\n", 
argv[0]);
exit(1);
}	
n = 1<<atoi(argv[1]);
m = 1<<atoi(argv[2]);
k = 1<<atoi(argv[3]); 
SIZE1= 1<<atoi(argv[4]);
SIZE2= 1<<atoi(argv[5]);


numTasks=1;
srand((unsigned) time(&timee)+200*rank);  
float tmp;
int i,cc,cout,rc;
float tx,ty,tz;  //Grid's boxes coordinates
int size=(SIZE1+SIZE2)-1;
int	pointNum=(SIZE1+SIZE2);

point *arr=(point *)malloc((2*pointNum)*sizeof(point));
int pointsPerBox=pointNum/(n*m*k)+32;

//set number 1	alliws... 1 array,size(2n).....1o miso->Q,2o ->C
point *a[n][m][k];
int counter[n][m][k];
for(i=0;i<n;i++){
  for(j=0;j<m;j++){
   for(cc=0;cc<k;cc++){
	a[i][j][cc]=(point *)malloc(pointsPerBox*sizeof(point));
	counter[i][j][cc]=0;
}
}
}
	
int sum=0;

	gettimeofday(&start0, NULL);
	 

for(i=0;i<pointNum;i++){
		sum+=1;
		tmp=((float)rand())/(RAND_MAX);
		arr[i].x=tmp;
		tmp=((float)rand())/(RAND_MAX);
		arr[i].y=tmp;
		tmp=((float)rand())/(RAND_MAX);
		arr[i].z=tmp;

		//	printf("PROCESS 0: x,y,z = %f,%f,%f \n",arr[i].x,arr[i].y,arr[i].z);	
		//	printf("PROCESS 1: x,y,z = %f,%f,%f \n",arr[i].x,arr[i].y,arr[i].z);
		if(i<=(pointNum/2))
			arr[i].type=1;
		else
			arr[i].type=2;
		//printf("%.2f\n",tmp);
		}
		//miinima....
		//checking and positioning setB in grid's boxes..
		// arr[i].process=rank;  
		float n1=(float)n;
		float m1=(float)m;
		float k1=(float)k;

		for(i=0;i<pointNum;i++){  //start for
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
		int t2=(int)ty;//printf("%d\n",t1); //point belongs
		int t3=(int)tz;  //increasing counter(frequency..)
		//	printf("lala %d\n",t1);
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
		if(rank==0)
		//	printf(" bz show....%d\n",((a[t1][t2][t3])[cout]).bz); //testing
			
		((a[t1][t2][t3])[cout]).process=arr[i].process;
		//	counter[t1][t2][t3]+=1;	//increasing counter(frequency..)

}//end for
int max1=-2;
for(i=0;i<n;i++){
 for(j=0;j<m;j++){
 for(cc=0;cc<k;cc++){
	if(counter[i][j][cc]>max1)  
		max1=counter[i][j][cc];
}}}
/*create struct type for MPI */

int cout0=0,tag=221,tmpCout;  //tmpCout holds new number of cells for each grid's box
int cox,coy,coz,flag,newCout;

//*******************************************************************************new code*******************************************************************************

gettimeofday(&start0, NULL);

point *res;
float fdist;
int row=1;
int col=(3*(SIZE2)+2);

int *sendCout=malloc(numTasks*sizeof(int));

for(j=0;j<(pointNum);j++){
	if((arr[j].type=1)){	
	//First checking points that contains our local array
	measure(arr,a,counter,n,m,k,j,row,col,sendCout,numTasks);//Function for finding (LOCAL)  min distance(min distance return)
	
	
	}
}


//***********************start new code send-receive for measure


 //for(j=0;j<pointNum;j++){

//	if(arr[j].maxDist!=INIT_MINDIST)
//	printf("FOR RANK : %d , MIN.DISTANCE OF POINT %f %f %f  is %f  .NEIGHBOUR IS : %f %f %f \n",rank,arr[j].x,arr[j].y,arr[j].z,arr[j].maxDist,arr[j].nbx,arr[j].nby,arr[j].nbz	);
 //} 




//MPI_Type_free(&rowtype);
	gettimeofday(&end0, NULL);
	 printf("FOR FIRST KERNEL, VERSION WITHOUT SHARED MEMORY TIME IS : %ld microSeconds \n", ((end0.tv_sec * 1000000 + end0.tv_usec)
		  - (start0.tv_sec * 1000000 + start0.tv_usec)));


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
//Allocation of 2d array that stores points to be sent to other processes

int cx=arr[j].bx;
int cy=arr[j].by;				//&&&&&&&&INFINITE LOOP KAPOY!!!!!!! 
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
			numb=counter[px][py][pz];
			for(i=0;i<numb;i++){
			if((((a[px][py][pz])[i]).type==2)){
			dist=getDistance(((a[px][py][pz])[i]).x,((a[px][py][pz])[i]).y	,(((a[px][py][pz])[i]).z),arr[j]);
			//	printf("2nd Distance %f..\n",dist);
			if((dist<(arr[j].maxDist))&&(dist!=0)) {    // !! Ensure that its not checking point with its self (!)
				arr[j].maxDist=dist;
			}
			}
		}		
	
}

//This function finds neighbour processes and send them candidates for nearest neighbours (slices of points that have contact with those processes
//Also receives as described above


//**********************************************************************************************************
void sendToproc(int *forSend,int pntx,int pnty,int pntz,int position,int col,int row,int *sendCout){

	int num=sendCout[position];
//	if(rank==0)
//		printf("NUM=%d ---rank= %d--- position=%d  point : %d , %d , %d\n",sendCout[position],rank,position,pntx,pnty,pntz);
		
	float coordinates = pntx | (pnty << 10) | (pntz << 20);
	forSend[position*col+num]=coordinates;
	
/*	if(rank==0){
	 fprintf(stdout, "x = %d\n", (forSend[position*col+num] & 0x3FF));   fprintf(stdout, "y = %d\n", ((forSend[position*col+num] >> 10) & 0x3FF));    fprintf(stdout, "z = %d\n", ((forSend[position*col+num] >> 20) & 0x3FF));
	}*/

	sendCout[position]+=1;
}
//**********************************************************************************************************


