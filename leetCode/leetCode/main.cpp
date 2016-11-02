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
	double pp[hidenode];//隐含结点的校正误差  
    double qq[outnode];//希望输出值与实际输出值的偏差  
    double yd[outnode];//希望输出值  
  
    double x[innode]; //输入向量  
    double x1[hidenode];//隐含结点状态值  
    double x2[outnode];//输出结点状态值  
    double o1[hidenode];//隐含层激活值  
    double o2[hidenode];//输出层激活值  
  
    for(int isamp=0;isamp<trainsample;isamp++)//循环训练一次样品  
    {  
        for(int i=0;i<innode;i++)  
            x[i]=p[isamp][i]; //输入的样本  
        for(int i=0;i<outnode;i++)  
            yd[i]=t[isamp][i]; //期望输出的样本  
  
        //构造每个样品的输入和输出标准  
        for(int j=0;j<hidenode;j++)  
        {  
            o1[j]=0.0;  
            for(int i=0;i<innode;i++)  
                o1[j]=o1[j]+weight1[i][j]*x[i];//隐含层各单元输入激活值  
            x1[j]=1.0/(1+exp(-o1[j] - threshold1[j]));//隐含层各单元的输出  
			//x1[j] = (o1[j] + threshold1[j]) >0 ?1:0;
        }  
        for(int k=0;k<outnode;k++)  
        {  
            o2[k]=0.0;  
            for(int j=0;j<hidenode;j++)  
                o2[k]=o2[k]+weight2[j][k]*x1[j]; //输出层各单元输入激活值  
            x2[k]=1.0/(1.0+exp(-o2[k] - threshold2[k])); //输出层各单元输出 
			//x2[k] = (o2[k]+ threshold2[k])>0?1:0;
		}  
  
        for(int k=0;k<outnode;k++)  
        {  
            qq[k]=(yd[k]-x2[k])*x2[k]*(1-x2[k]); //希望输出与实际输出的偏差  
            for(int j=0;j<hidenode;j++)  
                weight2[j][k]+=rate_w2*qq[k]*x1[j];  //下一次的隐含层和输出层之间的新连接权  
        }  
  
        for(int j=0;j<hidenode;j++)  
        {  
            pp[j]=0.0;  
            for(int k=0;k<outnode;k++)  
                pp[j]=pp[j]+qq[k]*weight2[j][k];  
            pp[j]=pp[j]*x1[j]*(1-x1[j]); //隐含层的校正误差  
  
            for(int i=0;i<innode;i++)  
                weight1[i][j]+=rate_w1*pp[j]*x[i]; //下一次的输入层和隐含层之间的新连接权  
        }  
        for(int k=0;k<outnode;k++)  
            e+=(yd[k]-x2[k])*(yd[k]-x2[k]); //计算均方差
        error=e/2.0;
        for(int k=0;k<outnode;k++)  
            threshold2[k]=threshold2[k]+rate_t2*qq[k]; //下一次的隐含层和输出层之间的新阈值  
        for(int j=0;j<hidenode;j++)  
            threshold1[j]=threshold1[j]+rate_t1*pp[j]; //下一次的输入层和隐含层之间的新阈值  
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
            active1[i]=active1[i]+weight1[j][i]*input[j]; //隐含层各单元激活值  
		hide1[i]=1.0/(1.0+exp(-active1[i]-threshold1[i])); //隐含层各单元输出 
	}
	for(int i=0;i<outnode;i++)
	{
		active2[i]=0.0;  
        for(int j=0;j<innode;j++)  
            active2[i]=active2[i]+weight2[j][i]*hide1[j]; //隐含层各单元激活值  
		hide2[i]=1.0/(1.0+exp(-active2[i]-threshold2[i])); //隐含层各单元输出 
	}
	for(int k=0;k<outnode;k++)
        result[k]=hide2[k];
	return result; 
}

//输入样本  
double X[trainsample][innode]= {  
    {0,0,0},{0,0,1},{0,1,0},{0,1,1},{1,0,0},{1,0,1},{1,1,0},{1,1,1}  
    };  
//期望输出样本  
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
	printf ("Please enter the no. of iterations\n请输入您要设定的迭代数 : ");  
    cin>>num; // 输入迭代次数 
    srand(time(0));    
    evpop(popcurrent );
    Max = popcurrent[0].fit;
    for(i =0;i< num;i ++)
    {
		cout<<"\ntimes ="<<i<<endl;
        for(j =0;j<4; j++)
            popnext[j ]=popcurrent[ j];  
        pickchrom(popnext ); // 挑选
        crossover(popnext ); // 交叉  
        mutation(popnext );  // 变异  
        for(j =0;j<4; j++)
            popcurrent[j ]=popnext[ j];// 种群更替；
    }  
	// 等待迭代终止；  
	//对于真正随机数是需要注意取较大的迭代次数  
    for(l =0;l<3; l++)  
    {  
        if(popcurrent [l]. fit > Max )  
        {  
            Max=popcurrent [l]. fit;  
            k=x(popcurrent [l]);//此时的value即为所求的x值  
        }
    }
    printf("\n 当x等于 %d时，函数得到最大值为： %d ",k ,Max);  
    printf("\nPress any key to end ! " );
    flushall();                                 // 清除所有缓冲区；  
    getche();                                   // 从控制台取字符，不以回车为结束；
	return 0;
}