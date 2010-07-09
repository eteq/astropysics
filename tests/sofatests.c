#include "sofa.h"
#include "stdio.h"
#include "stdlib.h"
#include "strings.h"
#define PI 3.1415926535897931

void printmat(char* title, double mat[3][3]) {
    printf("%s:\n[[ %e,%e,%e]\n [%e,%e,%e]\n[%e,%e,%e]]\n",title,
            mat[0][0],mat[0][1],mat[0][2],
            mat[1][0],mat[1][1],mat[1][2],
            mat[2][0],mat[2][1],mat[2][2]);
}

void printvec(char* title, double vec[3]) {
    printf("%s:\n[%e,%e,%e]\n",title,
            vec[0],vec[1],vec[2]);
}

void trans_matricies(double ep, double jd1) {
    double dpsi,deps,epsa,rb[3][3],rp[3][3],rbp[3][3],rn[3][3],rbpn[3][3],rc2i[3][3];

    iauNut00b(jd1,0,&dpsi,&deps);
    iauPn00(jd1,0,dpsi,deps,&epsa,rb,rp,rbp,rn,rbpn);
    printmat("B",rb);
    printmat("P",rp);
    printmat("N",rn);
    printmat("BPN",rbpn);

    iauC2i00b(jd1,0,rc2i);
    printmat("C",rc2i);

}

void trans_coords(double ep, double jd) {
    double vec[3],transvec[3],rc2i[3][3];
    double theta,phi,r;
    
    double lat = 10;
    double lng = 20;
    
    iauS2p(lng*PI/180,lat*PI/180,1,vec);
    iauC2i00b(jd,0,rc2i);
    printmat("C2I",rc2i);
    printvec("cartesian vec",vec);
    iauRxp(rc2i,vec,transvec);
    
    printvec("transvec",transvec);
    iauP2s(transvec,&phi,&theta,&r);
    
    double radeg = theta*180/PI/15;
    int rahr = (int)(radeg);
    int ramin = (int)((radeg-rahr)*60);
    double rasec = ((radeg-rahr-ramin/60.)*3600);
    double decdeg = phi*180/PI;
    int decdegi = (int)(decdeg);
    int decmin = (int)((decdeg-decdegi)*60);
    double decsec = ((decdeg-decdegi-decmin/60.)*3600);
    printf("transvec coords: ra:%ih%im%gs dec:%id%im%gs\n",rahr,ramin,rasec,decdegi,decmin,decsec);
    
    

}

int main(int argc, const char* argv[] ){
    if (argc != 4){
        printf("Need command line arguments: testname epoch jd\n");
    } else {
        double epoch,jd;
        char* testname = argv[1];
        int doall = strcmp(testname,"all")==0;
        
        epoch = atof(argv[2]);
        jd = atof(argv[3]);
        
        if (doall|(strcmp(testname,"trans_matricies")==0)) {
            trans_matricies(epoch,jd);
        }   
        if (doall|(strcmp(testname,"trans_coords")==0)) {
            trans_coords(epoch,jd);
        }  
    }
}
