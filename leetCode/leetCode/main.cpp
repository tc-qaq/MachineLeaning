#include <iostream>
#include <cmath>
#include <vector>
#include "approximate.h"

using namespace std;

#define innode 3
#define hidenode 10
#define outnode 1
#define trainsample 8
class BpNet
{
public :
	void train(double p[trainsample][innode],double t[trainsample][outnode]);
	double p[trainsample][innode];
	double t[trainsample][innode];
	double *recognize(double *p);
	BpNet();
    virtual ~BpNet();
public:
	void init();
	double weight1[innode][hidenode];
	double weight2[hidenode][outnode];
    double threshold1[hidenode];  
    double threshold2[outnode]; 
	double rate_w1;   
    double rate_w2;  
    double rate_t1;  
    double rate_t2;
	double e; 
    double error;  
    double result[outnode];  
};
BpNet::BpNet()
{
	error = 1.0;
	e = 0.0;
	rate_w1=0.9;  
    rate_w2=0.9;  
    rate_t1=0.9;   
    rate_t2=0.9; 
}
BpNet::~BpNet(){}
void winit(double w[],int n)
{
	for(int i=0;i<n;i++)
		w[i]=(2.0*(double)rand()/RAND_MAX)-1;
}
void BpNet::init()
{
	winit((double*)weight1,innode*hidenode);
	winit((double*)weight2,innode*hidenode);
	winit(threshold1,hidenode);
	winit(threshold2,hidenode);
}
void BpNet::train(double p[trainsample][innode],double t[trainsample][outnode])
{
	double pp[hidenode];//��������У�����  
    double qq[outnode];//ϣ�����ֵ��ʵ�����ֵ��ƫ��  
    double yd[outnode];//ϣ�����ֵ  
  
    double x[innode]; //��������  
    double x1[hidenode];//�������״ֵ̬  
    double x2[outnode];//������״ֵ̬  
    double o1[hidenode];//�����㼤��ֵ  
    double o2[hidenode];//����㼤��ֵ  
  
    for(int isamp=0;isamp<trainsample;isamp++)//ѭ��ѵ��һ����Ʒ  
    {  
        for(int i=0;i<innode;i++)  
            x[i]=p[isamp][i]; //���������  
        for(int i=0;i<outnode;i++)  
            yd[i]=t[isamp][i]; //�������������  
  
        //����ÿ����Ʒ������������׼  
        for(int j=0;j<hidenode;j++)  
        {  
            o1[j]=0.0;  
            for(int i=0;i<innode;i++)  
                o1[j]=o1[j]+weight1[i][j]*x[i];//���������Ԫ���뼤��ֵ  
            x1[j]=1.0/(1+exp(-o1[j] - threshold1[j]));//���������Ԫ�����  
			//x1[j] = (o1[j] + threshold1[j]) >0 ?1:0;
        }  
        for(int k=0;k<outnode;k++)  
        {  
            o2[k]=0.0;  
            for(int j=0;j<hidenode;j++)  
                o2[k]=o2[k]+weight2[j][k]*x1[j]; //��������Ԫ���뼤��ֵ  
            x2[k]=1.0/(1.0+exp(-o2[k] - threshold2[k])); //��������Ԫ��� 
			//x2[k] = (o2[k]+ threshold2[k])>0?1:0;
		}  
  
        for(int k=0;k<outnode;k++)  
        {  
            qq[k]=(yd[k]-x2[k])*x2[k]*(1-x2[k]); //ϣ�������ʵ�������ƫ��  
            for(int j=0;j<hidenode;j++)  
                weight2[j][k]+=rate_w2*qq[k]*x1[j];  //��һ�ε�������������֮���������Ȩ  
        }  
  
        for(int j=0;j<hidenode;j++)  
        {  
            pp[j]=0.0;  
            for(int k=0;k<outnode;k++)  
                pp[j]=pp[j]+qq[k]*weight2[j][k];  
            pp[j]=pp[j]*x1[j]*(1-x1[j]); //�������У�����  
  
            for(int i=0;i<innode;i++)  
                weight1[i][j]+=rate_w1*pp[j]*x[i]; //��һ�ε�������������֮���������Ȩ  
        }  
        for(int k=0;k<outnode;k++)  
            e+=(yd[k]-x2[k])*(yd[k]-x2[k]); //���������
        error=e/2.0;
        for(int k=0;k<outnode;k++)  
            threshold2[k]=threshold2[k]+rate_t2*qq[k]; //��һ�ε�������������֮�������ֵ  
        for(int j=0;j<hidenode;j++)  
            threshold1[j]=threshold1[j]+rate_t1*pp[j]; //��һ�ε�������������֮�������ֵ  
    } 
}
double *BpNet::recognize(double *p)
{
	double input[innode];
    double hide1[hidenode];
    double hide2[outnode];
    double active1[hidenode]; 
    double active2[hidenode];
	for(int i=0;i<innode;i++)
		input[i] = p[i];
	for(int i=0;i<hidenode;i++)
	{
		active1[i]=0.0;  
        for(int j=0;j<innode;j++)  
            active1[i]=active1[i]+weight1[j][i]*input[j]; //���������Ԫ����ֵ  
		hide1[i]=1.0/(1.0+exp(-active1[i]-threshold1[i])); //���������Ԫ��� 
	}
	for(int i=0;i<outnode;i++)
	{
		active2[i]=0.0;  
        for(int j=0;j<innode;j++)  
            active2[i]=active2[i]+weight2[j][i]*hide1[j]; //���������Ԫ����ֵ  
		hide2[i]=1.0/(1.0+exp(-active2[i]-threshold2[i])); //���������Ԫ��� 
	}
	for(int k=0;k<outnode;k++)
        result[k]=hide2[k];
	return result; 
}

//��������  
double X[trainsample][innode]= {  
    {0,0,0},{0,0,1},{0,1,0},{0,1,1},{1,0,0},{1,0,1},{1,1,0},{1,1,1}  
    };  
//�����������  
double Y[trainsample][outnode]={  
    {0},{0.125},{0.25},{0.375},{0.5},{0.625},{0.75},{0.875}  
    };  
int main()
{
	/*BpNet bp = BpNet();  
    bp.init();  
    int times=0;
	bp.error = 1;
	while(bp.error>0.0000000001)  
    {  
        bp.e=0.0;  
        times++;
        bp.train(X,Y);
		if(times%100000 == 0)
			cout<<"Times="<<times<<" error="<<bp.error<<endl;  
    }
	 cout<<"trainning complete..."<<endl;
	 double m[innode]={1,1,1};  
    double *r=bp.recognize(m);  
    for(int i=0;i<outnode;++i)  
       cout<<(bp.result[i])*8<<" ";
	system("pause");
	*/
	int num;
	int i,j,l,Max,k;
	Max = 0;
	printf ("Please enter the no. of iterations\n��������Ҫ�趨�ĵ����� : ");  
    cin>>num; // ����������� 
    srand(time(0));    
    evpop(popcurrent );
    Max = popcurrent[0].fit;
    for(i =0;i< num;i ++)
    {
		cout<<"\ntimes ="<<i<<endl;
        for(j =0;j<4; j++)
            popnext[j ]=popcurrent[ j];  
        pickchrom(popnext ); // ��ѡ
        crossover(popnext ); // ����  
        mutation(popnext );  // ����  
        for(j =0;j<4; j++)
            popcurrent[j ]=popnext[ j];// ��Ⱥ���棻
    }  
	// �ȴ�������ֹ��  
	//�����������������Ҫע��ȡ�ϴ�ĵ�������  
    for(l =0;l<3; l++)  
    {  
        if(popcurrent [l]. fit > Max )  
        {  
            Max=popcurrent [l]. fit;  
            k=x(popcurrent [l]);//��ʱ��value��Ϊ�����xֵ  
        }
    }
    printf("\n ��x���� %dʱ�������õ����ֵΪ�� %d ",k ,Max);  
    printf("\nPress any key to end ! " );
    flushall();                                 // ������л�������  
    getche();                                   // �ӿ���̨ȡ�ַ������Իس�Ϊ������
	return 0;
}