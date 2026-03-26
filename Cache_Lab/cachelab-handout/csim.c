#include "cachelab.h"
#include <getopt.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>

int clockk = 0 ;
int hits = 0 , misses = 0 , evictions = 0 ;
typedef struct Block
{
    int valid ;
    int time_stamp ;
    unsigned long tag ;
}Block;
int main(int argc , char* argv[])
{
    char c ; int s , E , b ; char* tracefile_path ;
    while((c = getopt(argc,argv,"s:E:b:t:"))!=-1) //!man 3 getopt来学习该函数用法
    {
        if(c=='s') s = atoi(optarg) ;
        else if(c=='E') E = atoi(optarg) ;
        else if(c=='b') b = atoi(optarg) ;
        else if(c=='t') tracefile_path = optarg ;
    }
    FILE* tracefile = fopen(tracefile_path,"r") ;
    int num_set = 1<<s ;
    Block* cache_memory = (Block *)malloc(num_set*E*sizeof(Block)) ;
    memset(cache_memory,0,sizeof(Block)*num_set*E) ;
    while(1)
    {
        unsigned long addr ; char op , isI ; 
        if(fscanf(tracefile,"%c%c %lx,%*d",&isI,&op,&addr)==EOF) break ;
        while(fscanf(tracefile,"%c",&c) && c!='\n') ;
        if(isI=='I') continue ;
        if(op=='M') hits++ ; //M,第二次一定hit
        unsigned long tag_bit = addr>>(s+b) , set_bit = (addr>>b)&((1<<s)-1) ;
        int idx = set_bit*E , lru = 0 ;
        for(int i=0 ; i<E ; i++)
        {
            Block* cnt_block = &cache_memory[idx+i] ;
            if(cnt_block->time_stamp<=cache_memory[idx+lru].time_stamp) lru = i ; //初始化timestamp==0,成立;逻辑上应该先找valid==0的block
            if(cnt_block->tag==tag_bit&&cnt_block->valid==1)
            {
                hits++ ; 
                cnt_block->time_stamp = ++clockk ;
                goto END ;
            }
        }
        //missed
        cache_memory[idx+lru].tag = tag_bit ; 
        cache_memory[idx+lru].time_stamp = ++clockk ;
        misses++ ;
        if(cache_memory[idx+lru].valid==1) evictions++ ;
        else cache_memory[idx+lru].valid = 1 ; 
        END:
    }
    printSummary(hits,misses,evictions); 
    fclose(tracefile) ;
    return 0;
}
