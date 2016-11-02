#include <iostream>
#include <stdio.h>  
#include <conio.h>  
#include <stdlib.h>  
#include <time.h>  
using namespace std;

typedef struct Chrom
{
	short int bit[6];
	int fit;
	double rate_fit;
	double conv_fit;
}chrom;

void evpop(chrom popcurrent[4]);//初始化
int x(chrom popcurrent); //转化
int y(int x);
void pickchrom(chrom popcurrent[4]);//选择
void pickchroms_new (chrom popcurrent[4]); // 基于概率分布
void crossover(chrom popcurrent[4]);//交叉
void mutation(chrom popcurrent[4]);//突变
double r8_uniform_ab(double a,double b,int &seed);//生成ab之间均匀的数字
chrom popcurrent[4];//初试种群
chrom popnext[4];//更新后种群

void evpop(chrom popcurrent[4])
{
	int i,j,value1;
	int random;
	double sum=0.0;
	for(i=0;i<4;i++)
	{
		for(j=0;j<6;j++)
		{
			random = rand();
			random %=2;
			popcurrent[i].bit[j] = random;
		}
		value1 = x(popcurrent[i]);
		popcurrent[i].fit = y(value1);
		sum += popcurrent[i].fit;
		// 输出整条染色体的编码情况
		printf("\n popcurrent[%d]=%d%d%d%d%d%d  value=%d  fitness = %d",i, popcurrent[i ].bit[5], popcurrent[i ].bit[4], popcurrent[i ].bit[3], popcurrent[i ].bit[2], popcurrent[i ].bit[1], popcurrent[i ].bit[0], value1,popcurrent [i]. fit);   
 	}
	for(i = 0;i<4;i++)
	{
		popcurrent[i].rate_fit = popcurrent[i].fit/sum;
		popcurrent[i].conv_fit = 0;
	}
}

int x(chrom popcurrent)
{
	int z;
	z = (popcurrent .bit[0]*1)+( popcurrent.bit [1]*2)+(popcurrent. bit[2]*4)+(popcurrent .bit[3]*8)+( popcurrent.bit [4]*16);
	if(popcurrent .bit[5]==1)  // 考虑到符号；  
    {  
        z=z *(-1);                               
    }
	return z;
}

int y(int x)
{
	int y;
	y = 5-x*x;
	return y;
}

void pickchrom(chrom popnext[4])
{
	int i,j;
	chrom temp;
	// 根据个体适应度来排序；（冒泡法）
	for(i =0;i<3; i++)  
    {  
        for(j =0;j<3-i; j++)  
        {  
            if(popnext [j+1]. fit>popnext [j]. fit)  
            {  
                temp=popnext [j+1];  
                popnext[j +1]=popnext[ j];  
                popnext[j ]=temp;
            }    
        }                 
    }
	for(i =0;i<4; i++)  
    {  
        printf("\nSorting:popnext[%d] fitness=%d" ,i, popnext[i ].fit);  
        printf("\n" );                       
    }  
    flushall();
}

void pickchroms_new(chrom popnext[4])
{
	int men,i,j;
	double p,sum=0.0;
	for (men = 0; men < 4; men++ )
        sum = sum + popnext[men].fit;//累积当前匹配值
	for (men = 0; men < 4; men++ )
        popnext[men].rate_fit = popnext[men].fit / sum;//计算每一个匹配率
	popcurrent[0].conv_fit = popcurrent[0].rate_fit;  //计算累积概率
    for ( men = 1; men < 4; men++)
        popnext[men].conv_fit = popnext[men-1].conv_fit + popnext[men].rate_fit;
	for ( i = 0; i < 4; i++ )  
    {
		//产生0~1之间的随机数  
        //p = r8_uniform_ab ( 0, 1, seed );//通过函数生成0~1之间均匀分布的数字  
        p =rand()%10;//  
        p = p/10;  
        if ( p < popnext[0].conv_fit )
            popcurrent[i] = popnext[0];
        else
            for ( j = 0; j < 4; j++ )
                if ( popnext[j].conv_fit <= p && p < popnext[j+1].conv_fit )
                    popcurrent[i] = popcurrent[j+1]; 
    }
	for ( i = 0; i < 4; i++ )
        popnext[i] = popcurrent[i];
}

double r8_uniform_ab(double a,double b,int &seed)
{
	int i4_huge = 2147483647;  
    int k;  
    double value;
	if ( seed == 0 )  
    {  
        std::cerr << "\n"; 
		std::cerr << "R8_UNIFORM_AB - Fatal error!\n";  
        std::cerr << "  Input value of SEED = 0.\n";  
        exit ( 1 );  
    }
	k = seed / 127773;  
    seed = 16807 * ( seed - k * 127773 ) - k * 2836;  
    if ( seed < 0 )
        seed = seed + i4_huge;
    value = ( double ) ( seed ) * 4.656612875E-10;  
    value = a + ( b - a ) * value;  
    return value;
}

void crossover(chrom popnext[4])
{
	int random, i;  
    random=(((int)rand() %5)+1);// 随机交叉点控制在0到5之间；  
    //0-rand (0,2)(1,3)
	for(i =0;i< random;i ++)         
    {
        popnext[2].bit [i]= popnext[0].bit [i];   // child 1 cross over  
        popnext[3].bit [i]= popnext[1].bit [i];   // child 2 cross over  
    }  
	//random-5 (0,3)(1,2)
    for(i =random; i<6;i ++) 
    {  
        popnext[2].bit [i]= popnext[1].bit [i];// child 1 cross over  
        popnext[3].bit [i]= popnext[0].bit [i];// chlid 2 cross over  
    }    
	// 为新个体计算适应度值
    for(i =0;i<4; i++)
        popnext[i ].fit= y(x (popnext[ i]));          
    // 输出新个体
	for(i =0;i<4; i++)
        printf("\nCrossOver popnext[%d]=%d%d%d%d%d%d    value=%d    fitness = %d",i, popnext[i ].bit[5], popnext[i ].bit[4], popnext[i ].bit[3], popnext[i ].bit[2], popnext[i ].bit[1], popnext[i ].bit[0], x(popnext [i]), popnext[i ].fit);   
}

void mutation(chrom popnext[4])
{
	int random,row ,col, value;
    random=rand ()%50;//0-50之间的数
    if(random ==25)  // 概率2%  
    {  
        col=rand ()%6; // 随机产生基因位号  
        row=rand ()%4; // 随机产生染色体号
		popnext[row].bit[col] = popnext[row].bit[col] == 0?1:0;//变异
		value=x (popnext[ row]);
        popnext[row].fit= y(value);
        printf("\nMutation occured in popnext[%d] bit[%d]:=%d%d%d%d%d%d    value=%d   fitness=%d", row,col ,popnext[ row].bit [5],popnext[ row].bit [4],popnext[ row].bit [3],popnext[ row].bit [2],popnext[ row].bit [1],popnext[ row].bit [0],value, popnext[row ].fit);  
    }
}
