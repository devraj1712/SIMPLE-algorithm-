#include<stdio.h>
#include<stdlib.h>
#include<math.h>

void initialize_vars(int M,int N,double u_prev[M+1][N][1],double v_prev[M][N+1][1], double u_star[M+1][N][1],double v_star[M][N+1][1],double p_star[M][N][1], double stream[M][N][1]);
double max_val(double a, double b, double c);
void solve_u_momentum(int M,int N,double u_prev[M+1][N][1],double v_prev[M][N+1][1], double u_star[M+1][N][1],double d[2*M-1][2*N-1][1],double p_star[M][N][1],int Rey_no);
void solve_v_momentum(int M,int N,double u_prev[M+1][N][1],double v_prev[M][N+1][1], double v_star[M][N+1][1],double d[2*M-1][2*N-1][1],double p_star[M][N][1],int Rey_no);
void solve_p_correctn(int M,int N,double u_star[M+1][N][1],double v_star[M][N+1][1], double p_dash[M][N][1], double d[2*M-1][2*N-1][1]);
int update_underrelax(int M, int N, int iter, double u_star[M+1][N][1],double v_star[M][N+1][1],double p_star[M][N][1],double p_dash[M][N][1],double u[M+1][N][1],double v[M][N+1][1],double d[2*M-1][2*N-1][1], double inf_norm[2], double u_prev[M+1][N][1],double v_prev[M][N+1][1]);
void compute_stream_fun(int M, int N,double stream[M][N][1], double u_star[M+1][N][1]);
void compute_vorticity(int M, int N,double u_star[M+1][N][1], double v_star[M][N+1][1], double vortic[M][N][1]);
void print_file(int M, int N, double stream[M][N][1], double u_star[M+1][N][1], double v_star[M][N+1][1], double vortic[M][N][1]);

int main()
{
    int M,N,Rey_no,i,j,iteration;
    printf("Enter the value of grid size MxN\n M = ");
    scanf("%d",&M);
    printf(" N = ");
    scanf("%d",&N);
    double u_star[M+1][N][1],v_star[M][N+1][1],u_prev[M+1][N][1],v_prev[M][N+1][1],d[2*M-1][2*N-1][1],p_star[M][N][1],p_dash[M][N][1];
    double u[M+1][N][1],v[M][N+1][1],inf_norm[2],stream[M][N][1],vortic[M][N][1];
    printf("Enter the value of Reynold's Number \n");
    scanf("%d",&Rey_no);
    iteration= 0;
    //initializing variables p,u,v
    initialize_vars(M, N, u_prev, v_prev, u_star, v_star, p_star, stream);
    do
    {
        if(iteration>0)
        {
            for(i=0;i<M;i++)
            {
                for(j=0;j<N;j++)
                {
                    u_prev[i][j][1] = u_star[i][j][1];
                    v_prev[i][j][1] = v_star[i][j][1];
                }
            }
        }
        //Solving u momentum equation
        solve_u_momentum(M,N,u_prev,v_prev,u_star,d,p_star,Rey_no);
        //Solving v momentum equation
        solve_v_momentum(M,N,u_prev,v_prev,v_star,d,p_star,Rey_no);
        //finding pressure correction terms
        solve_p_correctn(M,N,u_star,v_star,p_dash,d);
        //Using under relaxation to update the variables
        iteration = update_underrelax(M, N, iteration, u_star, v_star, p_star, p_dash, u, v, d, inf_norm, u_prev,v_prev);
    }while(inf_norm[0] > 0.000001 || inf_norm[1] > 0.000001);
    printf("\nSolution converged !");
    compute_stream_fun(M, N,stream, u_star);
    compute_vorticity(M, N, u_star, v_star, vortic);
    print_file(M, N, stream, u_star, v_star, vortic);
    return 0;
}

void initialize_vars(int M,int N,double u_prev[M+1][N][1],double v_prev[M][N+1][1], double u_star[M+1][N][1],double v_star[M][N+1][1],double p_star[M][N][1], double stream[M][N][1])
{
    for(int i=0;i<=M;i++)
    {
        for(int j=0;j<=N;j++)
        {
            if(i<M && j<N)
            {
                p_star[i][j][1] = 0.0;
                stream[i][j][1] = 0.0;
            }
            if(i<=M && j<N)
            {
                u_prev[i][j][1] =  0.0;
                u_star[i][j][1] =  0.0;
                //printf("%lf  ",u_prev[i][j]);
            }
            if(i<M && j<=N)
            {
                v_prev[i][j][1] = 0.0;
                v_star[i][j][1] = 0.0;
            }
        }
       // printf("\n");
    }
}

double max_val(double a, double b, double c)
{
    if (a>b)
    {
        if(a>c)
        {
            return a;
        }
        else
        {
            return c;
        }
    }
    else
    {
        if(b>c)
        {
            return b;
        }
        else
        {
            return c;
        }
    }
}

void solve_u_momentum(int M,int N,double u_prev[M+1][N][1],double v_prev[M][N+1][1], double u_star[M+1][N][1],double d[2*M-1][2*N-1][1], double p_star[M][N][1],int Rey_no)
{
    int i,j;
    double aP,aN,aS,aE,aW,Fn,Fs,Fe,Fw,Dn,Ds,De,Dw,del_x,del_y,Sp,Su;
    del_x = 1.0/(M);
    del_y = 1.0/(N);
    De = Dw = 1.0/(Rey_no*del_x);
    Dn = Ds = 1.0/(Rey_no*del_y);
    for(i=1;i<=M-1;i++)
    {
        for(j=0;j<=N-1;j++)
        {
            Fw = (u_prev[i-1][j][1] + u_prev[i][j][1])/2;
            Fe = (u_prev[i+1][j][1] + u_prev[i][j][1])/2;
            Fn = (v_prev[i-1][j+1][1] + v_prev[i][j+1][1])/2;
            Fs = (v_prev[i-1][j][1] + v_prev[i][j][1])/2;
            aW = max_val(Fw,Dw+Fw/2,0);
            aE = max_val(-1*Fe,De-Fe/2,0);
            aN = max_val(-1*Fn,Dn-Fn/2,0);
            aS = max_val(Fs,Ds+Fs/2,0);
            if(j==N-1)
            {
                Sp = -1.0/(Rey_no*del_y);
                Su = 1.0/(Rey_no*del_y);
                aP = aW + aE + aN + aS + Fe + Fn - Fw - Fs - Sp;
            }
            else if(j==0)
            {
                Sp = -1.0/(Rey_no*del_y);
                Su = 0;
                aP = aW + aE + aN + aS + Fe + Fn - Fw - Fs - Sp;
            }
            else
            {
                Su = 0;
                aP = aW + aE + aN + aS + Fe + Fn - Fw - Fs;
            }
            d[2*i-1][2*j][1] = 1.0/aP;
            u_star[i][j][1] = (aE*u_star[i+1][j][1] + aW*u_star[i-1][j][1] + aN*u_star[i][j+1][1] + aS*u_star[i][j-1][1] + p_star[i-1][j][1] - p_star[i][j][1] + Su)/aP;

            //printf("i=%d j=%d d=%lf\n",i,j,d[2*i-1][2*j][1]);
        }
    }
}

void solve_v_momentum(int M,int N,double u_prev[M+1][N][1],double v_prev[M][N+1][1], double v_star[M][N+1][1],double d[2*M-1][2*N-1][1], double p_star[M][N][1],int Rey_no)
{
    int i,j;
    double aP,aN,aS,aE,aW,Fn,Fs,Fe,Fw,Dn,Ds,De,Dw,del_x,del_y,Sp,Su;
    del_x = 1.0/(M);
    del_y = 1.0/(N);
    De = Dw = 1.0/(Rey_no*del_x);
    Dn = Ds = 1.0/(Rey_no*del_y);
    for(i=0;i<=M-1;i++)
    {
        for(j=1;j<=N-1;j++)
        {
            Fw = (u_prev[i][j][1] + u_prev[i][j-1][1])/2;
            Fe = (u_prev[i+1][j][1] + u_prev[i+1][j-1][1])/2;
            Fn = (v_prev[i][j+1][1] + v_prev[i][j][1])/2;
            Fs = (v_prev[i][j][1] + v_prev[i][j-1][1])/2;
            aW = max_val(Fw,Dw+Fw/2,0);
            aE = max_val(-1*Fe,De-Fe/2,0);
            aN = max_val(-1*Fn,Dn-Fn/2,0);
            aS = max_val(Fs,Ds+Fs/2,0);
            if(i==0 || i==M-1)
            {
                Sp = -1.0/(Rey_no*del_x);
                aP = aW + aE + aN + aS + Fe + Fn - Fw - Fs - Sp;
            }
            else
            {
                aP = aW + aE + aN + aS + Fe + Fn - Fw - Fs;
            }
            d[2*i][2*j-1][1] = 1.0/aP;
            v_star[i][j][1] = (aE*v_star[i+1][j][1] + aW*v_star[i-1][j][1] + aN*v_star[i][j+1][1] + aS*v_star[i][j-1][1] + p_star[i][j-1][1] - p_star[i][j][1] )/aP;

        }
    }
    return 0;
}

void solve_p_correctn(int M,int N,double u_star[M+1][N][1],double v_star[M][N+1][1],double p_dash[M][N][1],double d[2*M-1][2*N-1][1])
{
    int i,j;
    double aP,aN,aS,aE,aW;
    for(i=0;i<M;i++)
    {
        for(j=0;j<N;j++)
        {
            p_dash[i][j][1] = 0.0;
        }
    }
    //Solving pressure correction equation
    for(i=0;i<M;i++)
    {
        for(j=0;j<N;j++)
        {
            aW = 0;
            aE = 0;
            aS = 0;
            aN = 0;
            if(i!=0)
            {
                aW = d[2*i-1][2*j][1];
            }
            if(i!=M-1)
            {
                aE = d[2*i+1][2*j][1];
            }
            if(j!=0)
            {
                aS = d[2*i][2*j-1][1];
            }
            if(j!=N-1)
            {
                aN = d[2*i][2*j+1][1];
            }
            p_dash[i][j][1] = (aE*p_dash[i+1][j][1] + aW*p_dash[i-1][j][1] + aN*p_dash[i][j+1][1] + aS*p_dash[i][j-1][1] + u_star[i][j][1] - u_star[i+1][j][1] + v_star[i][j][1] - v_star[i][j+1][1])/(aE + aW + aN + aS);

        }
    }
}

int update_underrelax(int M, int N, int iter, double u_star[M+1][N][1],double v_star[M][N+1][1],double p_star[M][N][1],double p_dash[M][N][1],double u[M+1][N][1],double v[M][N+1][1],double d[2*M-1][2*N-1][1], double inf_norm[2],double u_prev[M+1][N][1],double v_prev[M][N+1][1])
{
    double alpha;
    inf_norm[0] = 0.0;
    inf_norm[1] = 0.0;
    alpha = 0.5;
    for(int i=0;i<M;i++)
    {
        for(int j=0;j<N;j++)
        {
            p_star[i][j][1] = p_star[i][j][1] + alpha*p_dash[i][j][1];
            if(i<M-1)
            {
                if(iter>0)
                {
                    u_prev[i+1][j][1] = u[i+1][j][1];
                }
                /*else
                {
                    u_prev[i+1][j][1] = 0.0;
                }*/
                u[i+1][j][1] = u_star[i+1][j][1] + d[2*i+1][2*j][1]*(p_dash[i][j][1] - p_dash[i+1][j][1]);
                if(fabs(u_prev[i+1][j][1] - u[i+1][j][1]) > inf_norm[0])
                {
                    inf_norm[0] = fabs(u_prev[i+1][j][1] - u[i+1][j][1]);
                }
                u_star[i+1][j][1] = alpha*u[i+1][j][1] + (1 - alpha)*u_prev[i+1][j][1];
            }
            if(j<N-1)
            {
                if(iter>0)
                {
                    v_prev[i][j+1][1] = v[i][j+1][1];
                }
                /*else
                {
                    v_prev[i][j+1][1] = 0.0;
                }*/
                v[i][j+1][1] = v_star[i][j+1][1] + d[2*i][2*j+1][1]*(p_dash[i][j][1] - p_dash[i][j+1][1]);
                if(fabs(v_prev[i][j+1][1] - v[i][j+1][1]) > inf_norm[1])
                {
                    inf_norm[1] = fabs(v_prev[i][j+1][1] - v[i][j+1][1]);
                }
                v_star[i][j+1][1] = alpha*v[i][j+1][1] + (1 - alpha)*v_prev[i][j+1][1];
            }
        }
    }
    iter++;
    printf("iteration %d u_l_inf %lf v_l_inf %lf\n",iter,inf_norm[0],inf_norm[1]);
    return iter;
}

void compute_stream_fun(int M, int N,double stream[M][N][1], double u_star[M+1][N][1])
{
    double norm_psi, old_value,del_y;
    del_y = 1.0/(N);
    do
    {
        norm_psi=0.0;
        for(int j=1;j<N;j++)
        {
            for(int i=1;i<M;i++)
            {
                old_value = stream[i][j][1];
                stream[i][j][1] = stream[i][j-1][1] + (u_star[i][j][1] + u_star[i+1][j][1])/2*del_y;
                if(fabs(stream[i][j][1] - old_value) > norm_psi)
                {
                    norm_psi = fabs(stream[i][j][1] - old_value);
                }
            }
        }
        //printf("\n norm_psi %lf",norm_psi);
    }while(norm_psi > .00001);
}

void compute_vorticity(int M, int N, double u_star[M+1][N][1], double v_star[M][N+1][1], double vortic[M][N][1])
{
    double del_x,del_y;
    del_x = 1.0/(M);
    del_y = 1.0/(N);
    for(int j=0;j<N;j++)
    {
        for(int i=0;i<M;i++)
        {
            //stream[i][j][1] = ( (u_star[i][j][1] + u_star[i+1][j][1])*del_y - (v_star[i][j][1] + v_star[i][j+1][1])*del_x )/2.0;
            if(i==0)
            {
                if(j==0)
                {
                    vortic[i][j][1] = ( (v_star[i][j+1][1] + v_star[i][j][1] + v_star[i+1][j+1][1] + v_star[i+1][j][1])/del_x - (u_star[i][j+1][1] + u_star[i+1][j+1][1] + u_star[i][j][1] + u_star[i+1][j][1])/del_y )/4.0;
                }
                else if(j==N-1)
                {
                    vortic[i][j][1] = ( (v_star[i][j+1][1] + v_star[i][j][1] + v_star[i+1][j+1][1] + v_star[i+1][j][1])/del_x - ( 4 - u_star[i][j][1] - u_star[i+1][j][1] - u_star[i][j-1][1] - u_star[i+1][j-1][1])/del_y )/4.0;
                }
                else
                {
                    vortic[i][j][1] = ( (v_star[i][j+1][1] + v_star[i][j][1] + v_star[i+1][j+1][1] + v_star[i+1][j][1])/del_x - (u_star[i][j+1][1] + u_star[i+1][j+1][1] - u_star[i][j-1][1] - u_star[i+1][j-1][1])/del_y )/4.0;
                }
            }
            else if(i==M-1)
            {
                if(j==0)
                {
                    vortic[i][j][1] = ( ( 0 - v_star[i-1][j+1][1] - v_star[i-1][j][1] - v_star[i][j+1][1] - v_star[i][j][1])/del_x - (u_star[i][j+1][1] + u_star[i+1][j+1][1] + u_star[i][j][1] + u_star[i+1][j][1])/del_y )/4.0;
                }
                else if(j==N-1)
                {
                    vortic[i][j][1] = ( ( 0 - v_star[i-1][j+1][1] - v_star[i-1][j][1] - v_star[i][j+1][1] - v_star[i][j][1])/del_x - ( 4 - u_star[i][j-1][1] - u_star[i+1][j-1][1] - u_star[i][j][1] - u_star[i+1][j][1])/del_y )/4.0;
                }
                else
                {
                    vortic[i][j][1] = ( ( 0 - v_star[i-1][j+1][1] - v_star[i-1][j][1] - v_star[i][j+1][1] - v_star[i][j][1])/del_x - (u_star[i][j+1][1] + u_star[i+1][j+1][1] - u_star[i][j-1][1] - u_star[i+1][j-1][1])/del_y )/4.0;
                }
            }
            else
            {
                if(j==0)
                {
                    vortic[i][j][1] = ( (v_star[i+1][j+1][1] + v_star[i+1][j][1] - v_star[i-1][j+1][1] - v_star[i-1][j][1])/del_x - (u_star[i][j+1][1] + u_star[i+1][j+1][1] + u_star[i][j][1] + u_star[i+1][j][1])/del_y )/4.0;
                }
                else if(j==N-1)
                {
                    vortic[i][j][1] = ( (v_star[i+1][j+1][1] + v_star[i+1][j][1] - v_star[i-1][j+1][1] - v_star[i-1][j][1])/del_x - ( 4 - u_star[i][j][1] - u_star[i+1][j][1] - u_star[i][j-1][1] - u_star[i+1][j-1][1])/del_y )/4.0;
                }
                else
                {
                    vortic[i][j][1] = ( (v_star[i+1][j+1][1] + v_star[i+1][j][1] - v_star[i-1][j+1][1] - v_star[i-1][j][1])/del_x - (u_star[i][j+1][1] + u_star[i+1][j+1][1] - u_star[i][j-1][1] - u_star[i+1][j-1][1])/del_y )/4.0;
                }
            }
        }
    }
}

void print_file(int M, int N, double stream[M][N][1], double u_star[M+1][N][1], double v_star[M][N+1][1], double vortic[M][N][1])
{
    FILE *f1=fopen("stream.plt", "w");
    FILE *f2=fopen("velocity.plt", "w");
    FILE *f3=fopen("u_profile.plt", "w");
    FILE *f4=fopen("v_profile.plt", "w");
    FILE *f5=fopen("vorticity.plt", "w");
    fprintf(f1,"VARIABLES =\"X\", \"Y\", \"PSI\"\n");
    fprintf(f1, "ZONE T = \"BLOCK1\", I = %d, J = %d, F = POINT\n\n",M,N);
    fprintf(f2,"VARIABLES =\"X\", \"Y\", \"U\", \"V\"\n");
    fprintf(f2, "ZONE T = \"BLOCK1\", I = %d, J = %d, F = POINT\n\n",M,N);
    fprintf(f3,"VARIABLES = \"Y\", \"U\"\n");
    fprintf(f3, "ZONE T = \"BLOCK1\", I = %d, F = POINT\n\n",N);
    fprintf(f4,"VARIABLES = \"X\", \"V\"\n");
    fprintf(f4, "ZONE T = \"BLOCK1\", I = %d, F = POINT\n\n",M);
    fprintf(f5,"VARIABLES =\"X\", \"Y\", \"Omega\"\n");
    fprintf(f5, "ZONE T = \"BLOCK1\", I = %d, J = %d, F = POINT\n\n",M,N);
    int i,j;
    double del_x,del_y;
    del_x = 1.0/(M);
    del_y = 1.0/(N);
    for(j=0;j<N;j++)
    {
        for(i=0;i<M;i++)
        {
            fprintf(f1, "%lf \t %lf \t %lf \n", i*del_x, j*del_y, stream[i][j][1]);
            fprintf(f2, "%lf \t %lf \t %lf \t %lf \n", i*del_x, j*del_y, u_star[i][j][1], v_star[i][j][1]);
            if(i==M/2)
            {
                fprintf(f3, " %lf \t %lf \n", j*del_y, u_star[i][j][1]);
            }
            if(j==N/2)
            {
                fprintf(f4, " %lf \t %lf \n", i*del_x, v_star[i][j][1]);
            }
            fprintf(f5, "%lf \t %lf \t %lf \n", i*del_x, j*del_y, vortic[i][j][1]);
        }
    }
    fclose(f1);
    fclose(f2);
    fclose(f3);
    fclose(f4);
    fclose(f5);
}
