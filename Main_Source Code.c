#include<stdio.h>
#include<math.h>
#include<conio.h>
#include<stdlib.h>
void trans(float [][6], float [][6], int,int,int,float *Ref3); //for transpose of matrix   A matrix is:
void cofac(float [][6], int,int,int,float*Ref3); //for cofactors of matrix
float determin(float [][6], int); //for determinant of matrix
int main()
{
 int i,j,k,a,an;
 float Ef1,Ef2,Gf,Em,Gm,Vf,vf,vm,E1,E2,G12,v12,Q[3][3];
 // Code for calculating effective properties started
 printf("Enter Different required values in Gpa\n");
 printf("Fiber longitudinal Modulus:");
 scanf("%f",&Ef1);
 printf("Fiber Transverse Modulus:");
 scanf("%f",&Ef2);
 printf("Matrix Modulus:");
 scanf("%f",&Em);
 printf("Fiber Shear Modulus:");
 scanf("%f",&Gf);
 printf("Matrix Shear Modulus:");
 scanf("%f",&Gm);
 printf("Fiber poisson ratio:");
 scanf("%f",&vf);
 printf("Matrix poisson ratio:");
 scanf("%f",&vm);
 printf("Fiber volume fraction:");
 scanf("%f",&Vf);
 E1=(Ef1*Vf)+(Em)*(1.0-Vf);
 E2=(Ef1*Em)/(Ef1*(1.0-Vf)+Em*Vf);
 G12=(Gf*Gm)/(Gf*(1.0-Vf)+Gm*Vf);
 v12=vf*Vf+vm*(1.0-Vf);
 printf("\nThe effective properties are:\n");
 printf("Longitudinal Modulus is:%0.3f Gpa\n",E1);
 printf("Transverse Modulus is:%0.3f Gpa2\n",E2);
 printf("Shear Modulus is:%0.3f Gpa\n",G12);
 printf("Poisson Ratio is:%0.3f\n",v12);
 // Code ended
 // Code for calculating Q matrix started
 printf("\nHow many laminates are present:");
 scanf("%d",&an);
 int Angle[an];
 printf("\nEnter 0,41,-41,90,90,-41,41,0 angles in degrees from top to bottom with pressing enter button at end again enter 1:\n");
 for(i=0;i<an;i++){
    scanf("%d",&Angle[i]);
 }
  for(i=0;i<=2;i++)
  {
    for(j=0;j<=2;j++)
    {
      if(i==0&&j==0)
      {
        Q[i][j]= (E1)/(1.0-v12*v12*(E2/E1));
      }
      else if(i==1&&j==1)
      {
        Q[i][j]= (E2)/(1.0-v12*v12*(E2/E1));
      }
      else if(i==0&&j==1||i==1&&j==0)
      {
        Q[i][j]=(v12*E2)/(1.0-v12*v12*(E2/E1));
      }
      else if(i==2&&j==2)
      {
        Q[i][j]=G12;
      }
      else
      {
        Q[i][j]=0;
      }
    }
  }
  printf("\nQ matrices for each angles:\n");
  float *Arr=(float*)malloc(9*an*sizeof(float));
  float *Ref=Arr;
  float *Ref3=Arr; // pointer to use Q values in stress calculation
  for(k=0;k<an;k++)
  {
    float angle= 22.0*Angle[k]/(7*180);
    float m=cos(angle),n=sin(angle);
    float *I=Arr;
    *Arr=pow(m,4)*Q[0][0]+2*pow(m,2)*pow(n,2)*(Q[0][1]+Q[2][2]*2)+pow(n,4)*Q[1][1];
    Arr++;
    *Arr=pow(m,2)*pow(n,2)*(Q[0][0]+Q[1][1]-4*Q[2][2])+ Q[0][1]*(pow(m,4)+pow(n,4));
    Arr++;
    *Arr=n*pow(m,3)*(Q[0][0]-Q[0][1])+m*pow(n,3)*(Q[0][1]-Q[1][1])-2*m*n*Q[2][2]*(m*m-n*n);
    Arr++;
    *Arr=pow(m,2)*pow(n,2)*(Q[0][0]+Q[1][1]-4*Q[2][2])+ Q[0][1]*(pow(m,4)+pow(n,4));
    Arr++;
    *Arr=pow(n,4)*Q[0][0]+2*pow(m,2)*pow(n,2)*(Q[0][1]+Q[2][2]*2)+pow(m,4)*Q[1][1];
    Arr++;
    *Arr=m*pow(n,3)*(Q[0][0]-Q[0][1])+n*pow(m,3)*(Q[0][1]-Q[1][1])+2*m*n*Q[2][2]*(m*m-n*n);
    Arr++;
    *Arr=n*pow(m,3)*(Q[0][0]-Q[0][1])+m*pow(n,3)*(Q[0][1]-Q[1][1])-2*m*n*Q[2][2]*(m*m-n*n);
    Arr++;
    *Arr=m*pow(n,3)*(Q[0][0]-Q[0][1])+n*pow(m,3)*(Q[0][1]-Q[1][1])+2*m*n*Q[2][2]*(m*m-n*n);
    Arr++;
    *Arr=pow(m,2)*pow(n,2)*(Q[0][0]+Q[1][1]-2*Q[0][1])+(m*m-n*n)*Q[2][2];
    Arr++;
    printf("\nQ matrix for angle %d:\n",Angle[k]);
    for(i=0;i<=2;i++){
        for(j=0;j<=2;j++){
            printf("%f ",*I);
            I++;
        }
        printf("\n");
    }
  }
  // Code Ended
  // Code to calculate A,B,D matrix Started
   printf("\nNow we will calculate A,B,D matrices:");
   float A[3][3],B[3][3],D[3][3];
   float T=1.0; // Total thickness in mm
   float t=(T*1.0)/an;// Thickness of each layer
   float *Ref2=Ref;
   for(i=0;i<3;i++){
    for(j=0;j<3;j++){
        A[i][j]= *Ref*(-3*t+4*t)+*(Ref+=9)*((-2+3)*t)+*(Ref+=9)*((-1+2)*t)+*(Ref+=9)*((0+1)*t)+*(Ref+=9)*((1-0)*t)+*(Ref+=9)*((2-1)*t)+*(Ref+=9)*((3-2)*t)+*(Ref+=9)*((4-3)*t);
        Ref=Ref2;
        B[i][j]= 0.5*(*Ref*((-3)*(-3)*t*t-(-4)*(-4)*t*t)+*(Ref+=9)*(((-2)*(-2)-(-3)*(-3))*t*t)+*(Ref+=9)*(((-1)*(-1)-(-2)*(-2))*t*t)+*(Ref+=9)*(((0)*(0)-(-1)*(-1))*t*t)+*(Ref+=9)*(((1)*(1)-(0)*(0))*t*t)+*(Ref+=9)*(((2)*(2)-(1)*(1))*t*t)+*(Ref+=9)*(((3)*(3)-(2)*(2))*t*t)+*(Ref+=9)*(((4)*(4)-(3)*(3))*t*t));
        Ref=Ref2;
        D[i][j]= (1.0/3)*(*Ref*(pow(-3*t,3)-pow(-4*t,3))+*(Ref+=9)*(pow(-2*t,3)-pow(-3*t,3))+*(Ref+=9)*(pow(-1*t,3)-pow(-2*t,3))+*(Ref+=9)*(pow(0*t,3)-pow(-1*t,3))+*(Ref+=9)*(pow(1*t,3)-pow(0*t,3))+*(Ref+=9)*(pow(2*t,3)-pow(1*t,3))+*(Ref+=9)*(pow(3*t,3)-pow(2*t,3))+*(Ref+=9)*(pow(4*t,3)-pow(3*t,3)));
        Ref=Ref2;
        Ref++;
        Ref2=Ref;
    }
   }
  // Adjusting -0.000   values to 0
  for(i=0;i<3;i++){
    for(j=0;j<3;j++){
        if(A[i][j]<1){
            A[i][j]=0;
        }
        if(B[i][j]<1){
            B[i][j]=0;
        }
        if(D[i][j]<1){
            D[i][j]=0;
        }
    }
  }
  printf("\nA matrix is:\n");
  for(i=0;i<=2;i++)
   {
    for(j=0;j<=2;j++)
     {
       printf("%0.3f  ",A[i][j]);
    }
    printf("\n");
   }
   printf("\nB matrix is:\n");
  for(i=0;i<=2;i++)
   {
    for(j=0;j<=2;j++)
     {
       printf("%0.3f  ",B[i][j]);
    }
    printf("\n");
   }
   printf("\nD matrix is:\n");
  for(i=0;i<=2;i++)
   {
    for(j=0;j<=2;j++)
     {
       printf("%0.3f  ",D[i][j]);
    }
    printf("\n");
   }
   // Code Ended
   // Calculating Inverse matrix of combined A,B,D matrix
    int  row=6, col=6;
    float matrix[6][6];
    for(i=0;i<6;i++){
        for(j=0;j<6;j++){
            if(i<3&&j<3){
                matrix[i][j]=A[i][j];
            }
            else if(i>=3&&j>=3){
                matrix[i][j]=D[i-3][j-3];
            }
            else{
                if(j<3){
                    matrix[i][j]=B[i-3][j];
                }
                else if(j>3){
                    matrix[i][j]=B[i][j-3];
                }
            }
        }
    }
   // Code ended
   // Calculating inverse matrix of combined matrix to get strain equation with applied forces and moments
   float inv_matrix[6][6];
   cofac(matrix, row,t,an,Ref3);
 return 0;
}

// This function is to find cofactors of matrix . . .
void cofac(float comatr[6][6], int f,int t,int an,float *Ref3){
    float matr[6][6], cofact[6][6];
    int a, b, c, d, x, y;
    for(c=0; c<f; ++c){
     for(d=0; d<f; ++d){
        x=0;
        y=0;
        for(a=0;a<f; ++a){
            for(b=0; b<f; ++b){
               if(a != c && b != d)
               {
               matr[x][y]=comatr[a][b];
               if(y<(f-2))
               y++;
               else
               {
               y=0;
               x++;
               }
               }
            }
        }
        cofact[c][d] = pow(-1,c + d) * determin(matr,f-1);
     }
    }
    trans(comatr, cofact, f,t,an,Ref3);
}
// This function is to find the determinant value of matrix . . .
float determin(float matrix[6][6], int k)
{
    float deter=0.0, z=1.0, mt[6][6];
    int a, b, c, x, y;
    if(k==1)
    {
      return(matrix[0][0]);
     }
    else
    {
        deter=0;
        for(a=0;a<k;++a){
        x=0;
        y=0;
        for(b=0;b<k;++b){
        for(c=0;c<k;++c){
        mt[b][c]=0;
        if((b != 0) && (c != a))
        {
         mt[x][y]=matrix[b][c];
         if(y<(k-2))
           y++;
         else
         {
            y=0;
            x++;
          }
         }
            }
            }
        deter=deter + z * (matrix[0][a] * determin(mt,k-1));
         z=-1 * z;
            }
        }
        return(deter);
}

// This function is to transpose of matrix . . .
void trans(float matr[6][6], float m1[6][6], int r,int t,int an,float *Ref3)
{
  float inv_matrix[6][6], tran[6][6], d;
  int a,b;
  for(a=0;a<r;++a){
    for(b=0;b<r;++b){
        tran[a][b]=m1[b][a];
    }
  }
  d=determin(matr, r);
  for(a=0;a<r;++a){
    for(b=0;b<r;++b){
        inv_matrix[a][b]=tran[a][b] / d;
    }
  }

     // Calculating laminate strains Strains
   float e_01,e_02,e_012; // in-plane strains
   float e1,e2,e12;       // Laminate strains
   float K1,K2,K12;       // Mid-Surface curvatures
   float N1,N2,N12,M1,M2,M12; // Forces and Moments
   printf("Enter Forces and moments in unit N/mm ,press enter after every value:\n");
   printf("Nx:");
   scanf("%f",&N1);
   printf("Ny:");
   scanf("%f",&N2);
   printf("Nxy:");
   scanf("%f",&N12);
   printf("Mx:");
   scanf("%f",&M1);
   printf("My:");
   scanf("%f",&M2);
   printf("Mxy:");
   scanf("%f",&M12);
    e_01=(inv_matrix[0][0]*N1+inv_matrix[0][1]*N2+inv_matrix[0][2]*N12+inv_matrix[0][3]*M1+inv_matrix[0][4]*M2+inv_matrix[0][5]*M12)*1.0;
    e_02=(inv_matrix[1][0]*N1+inv_matrix[1][1]*N2+inv_matrix[1][2]*N12+inv_matrix[1][3]*M1+inv_matrix[1][4]*M2+inv_matrix[1][5]*M12)*1.0;
    e_012=(inv_matrix[2][0]*N1+inv_matrix[2][1]*N2+inv_matrix[2][2]*N12+inv_matrix[2][3]*M1+inv_matrix[2][4]*M2+inv_matrix[2][5]*M12)*1.0;
    K1=(inv_matrix[3][0]*N1+inv_matrix[3][1]*N2+inv_matrix[3][2]*N12+inv_matrix[3][3]*M1+inv_matrix[3][4]*M2+inv_matrix[3][5]*M12)*1.0;
    K2=(inv_matrix[4][0]*N1+inv_matrix[4][1]*N2+inv_matrix[4][2]*N12+inv_matrix[4][3]*M1+inv_matrix[4][4]*M2+inv_matrix[4][5]*M12)*1.0;
    K12=(inv_matrix[5][0]*N1+inv_matrix[5][1]*N2+inv_matrix[5][2]*N12+inv_matrix[5][3]*M1+inv_matrix[5][4]*M2+inv_matrix[5][5]*M12)*1.0;
    // Laminate Strains
    e1=e_01+t*K1;
    e2=e_02+t*K2;
    e12=e_012+t*K12;
    // Calculating Stresses in each laminates
    float *Sigma_1=(float*)malloc(an*sizeof(float)); // Creating pointer to store longitudinal stresses in 8 laminates
    float *Sigma_2=(float*)malloc(an*sizeof(float)); // Creating pointer to store Transverse stresses in 8 laminates
    float *Sigma_12=(float*)malloc(an*sizeof(float)); // Creating pointer to store Shear stresses in 8 laminates
    int i,j;
    for(i=0;i<an;i++){
        *Sigma_1=(*Ref3)*e1+*(Ref3+=1)*e2+*(Ref3+=1)*e12;
        *Sigma_2=*(Ref3+=1)*e1+*(Ref3+=1)*e2+*(Ref3+=1)*e12;
        *Sigma_12=*(Ref3+=1)*e1+*(Ref3+=1)*e2+*(Ref3+=1)*e12;
        Ref3+=1;
        printf("\nStresses in lamina_%d:\n",i+1);
        printf("Sigma_x: %f\nSigma_y: %f\nSigma_xy: %f\n",*Sigma_1,*Sigma_2,*Sigma_12);
        Sigma_1++;
        Sigma_2++;
        Sigma_12++;
    }
    Sigma_1-=8;
    Sigma_2-=8;
    Sigma_12-=8;

    // Strength values in Gpa
    float X=1.447;
    float Xc=1.447;
    float Y=0.052;
    float Yc=0.206;
    float S_12=0.093;
    float S_23=0.034;
    // Converting Stresses from x-y axes to 1-2 axes
    // Writing angles again
    int Angle[8]={0,41,-41,90,90,-41,41,0};
    float Rad;
    float *sigma_1=(float*)malloc(an*sizeof(float)); // Creating pointer to store longitudinal stresses
    float *sigma_2=(float*)malloc(an*sizeof(float)); // Creating pointer to store Transverse stresses
    float *sigma_12=(float*)malloc(an*sizeof(float)); // Creating pointer to store Shear stresses
     printf("\nStresses in (1-2) Axis:\n");
    for(i=0;i<8;i++){
         Rad=22.0*Angle[i]/(7*180);
        *sigma_1=(*Sigma_1+*Sigma_2)/2.0+(*Sigma_1-*Sigma_2)*cos(2*Rad)/2.0+*Sigma_12*sin(2*Rad);
        *sigma_2=(*Sigma_1+*Sigma_2)/2.0-(*Sigma_1-*Sigma_2)*cos(2*Rad)/2.0-*Sigma_12*sin(2*Rad);
        *sigma_12=-(*Sigma_1-*Sigma_2)*sin(2*Rad)/2.0+*Sigma_12*cos(2*Rad);
        printf("\nStresses in Ply:%d",i+1);
        printf("\n   sigma_1: %0.3f Gpa,  sigma_2: %0.3f Gpa,  sigma_12: %0.3f Gpa",*sigma_1,*sigma_2,*sigma_12);
        Sigma_1++;
        Sigma_2++;
        Sigma_12++;
        sigma_1++;
        sigma_2++;
        sigma_12++;
    }
    sigma_1-=8; // Placing pointer to Starting address
    sigma_2-=8;
    sigma_12-=8;
    // Using Hashin failure criteria
    printf("\nUsing Hashin Failure Criteria:\n");
    float a1,b1,c1,r1,dis,R1;
    float C;
    for(i=0;i<an;i++){
        if(*sigma_2>0){
          C= (*sigma_2/Y)*(*sigma_2/Y)+(*sigma_12/S_12)*(*sigma_12/S_12);
          if(C>=1){
            printf("\nLamina %d fails by Matrix Tensile failure:\n",i+1);
            a1=(*sigma_2/Y)*(*sigma_2/Y)+(*sigma_12/S_12)*(*sigma_12/S_12);
            b1=0;
            c1=(-1);
            dis=(b1*b1)-(4.0*a1*c1);
            R1=((-b1)+sqrt(dis))/(2.0*a1);
            printf("Margin of Safety for lamina %d, R: %f",i+1,R1);
            sigma_1++;
            sigma_2++;
            sigma_12++;
            printf("\n");
            continue;

          }
        }
        else if (*sigma_2<0){
            C=(pow(*sigma_2/(4*S_23),2))+(pow(*sigma_12/S_12,2));
            if(C>=1){
            printf("\nLamina %d fails by Matrix Compressive failure:\n",i+1);
            a1=(pow(*sigma_2/(4*S_23),2))+(pow(*sigma_12/S_12,2));
            b1=(pow(Yc/(2.0*S_23),2)-1)*(*sigma_2/Yc);
            c1=(-1);
            dis=(b1*b1)-(4.0*a1*c1);
            R1=((-b1)+sqrt(dis))/(2.0*a1);
            printf("Margin of Safety for lamina %d, R: %f",i+1,R1);
            sigma_1++;
            sigma_2++;
            sigma_12++;
            printf("\n");
            continue;
          }
        }
        C=(*Sigma_1/X)*(*sigma_1/X)+(*sigma_12/S_12)*(*sigma_12/S_12);
        if(C>=1){
            printf("\nLamina %d fails by Tensile Fiber failure:\n",i+1);
            a1=(*Sigma_1/X)*(*sigma_1/X)+(*sigma_12/S_12)*(*sigma_12/S_12);
            b1=0;
            c1=(-1);
            dis=(b1*b1)-(4.0*a1*c1);
            R1=((-b1)+sqrt(dis))/(2.0*a1);
            printf("Margin of Safety for lamina %d, R: %f",i+1,R1);
            sigma_1++;
            sigma_2++;
            sigma_12++;
            printf("\n");
            continue;
        }
        C=(*sigma_1/(Xc*1.0))*(*sigma_1/(Xc*1.0));
        if(C>=1){
            printf("\nLamina %d fails by Compressive Fiber failure:\n",i+1);
            a1=(*sigma_1/(Xc*1.0))*(*sigma_1/(Xc*1.0));
            b1=0;
            c1=(-1);
            dis=(b1*b1)-(4.0*a1*c1);
            R1=((-b1)+sqrt(dis))/(2.0*a1);
             printf("Margin of Safety for lamina %d, R: %f",i+1,R1);
             sigma_1++;
             sigma_2++;
             sigma_12++;
             printf("\n");
             continue;
        }
        if(C<1){
            printf("\nLamina %d is safe\n",i+1);
            a1=(*sigma_2/Y)*(*sigma_2/Y)+(*sigma_12/S_12)*(*sigma_12/S_12);
            b1=0;
            c1=(-1);
            dis=(b1*b1)-(4.0*a1*c1);
            R1=((-b1)+sqrt(dis))/(2.0*a1);
            printf("Margin of Safety for lamina %d, R: %f",i+1,R1);
            sigma_1++;
            sigma_2++;
            sigma_12++;
            printf("\n");
            continue;
        }
    }
}
