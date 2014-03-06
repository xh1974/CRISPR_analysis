/*
 *  RRA.c
 *  Implementation of Robust Rank Aggregation (RRA)
 *
 *  Created by Han Xu on 20/12/13.
 *  Copyright 2013 Dana Farber Cancer Institute. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <assert.h>
#include "math_api.h"
#include "words.h"
#include "rvgs.h"
#include "rngs.h"

#define MAX_NAME_LEN 255           //maximum length of item name, group name or list name
#define CDF_MAX_ERROR 1E-10        //maximum error in Cumulative Distribution Function estimation in beta statistics
#define MAX_GROUP_NUM 100000       //maximum number of groups
#define MAX_LIST_NUM 1000          //maximum number of list 
#define RAND_PASS_NUM 100          //number of passes in random simulation for computing FDR

typedef struct
{
	char name[MAX_NAME_LEN];       //name of the item
	int listIndex;                 //index of list storing the item
	double value;                  //value of measurement
	double percentile;             //percentile in the list
} ITEM_STRUCT;

typedef struct 
{
	char name[MAX_NAME_LEN];       //name of the group
	ITEM_STRUCT *items;            //items in the group
	int itemNum;                   //number of items in the group
	double loValue;                //lo-value in RRA
	double fdr;                    //false discovery rate
} GROUP_STRUCT;

typedef struct
{
	char name[MAX_NAME_LEN];       //name of the list
	double *values;                //values of items in the list
	int itemNum;                   //number of items in the list
} LIST_STRUCT;

//Read input file. File Format: <item id> <group id> <list id> <value>. Return 1 if success, -1 if failure
int ReadFile(char *fileName, GROUP_STRUCT *groups, int maxGroupNum, int *groupNum, LIST_STRUCT *lists, int maxListNum, int *listNum);

//Save group information to output file. Format <group id> <number of items in the group> <lo-value> <false discovery rate>
int SaveGroupInfo(char *fileName, GROUP_STRUCT *groups, int groupNum);

//Process groups by computing percentiles for each item and lo-values for each group
int ProcessGroups(GROUP_STRUCT *groups, int groupNum, LIST_STRUCT *lists, int listNum, double maxPercentile);

//QuickSort groups by loValue
void QuickSortGroupByLoValue(GROUP_STRUCT *groups, int start, int end);

//Compute False Discovery Rate based on uniform distribution
int ComputeFDR(GROUP_STRUCT *groups, int groupNum, double maxPercentile, int numOfRandPass);

//print the usage of Command
void PrintCommandUsage(const char *command);


//Compute lo-value based on an array of percentiles
int ComputeLoValue(double *percentiles,     //array of percentiles
				   int num,                 //length of array
				   double *loValue,         //pointer to the output lo-value
				   double maxPercentile);   //maximum percentile, computation stops when maximum percentile is reached

int main (int argc, const char * argv[]) 
{
	int i,flag;
	GROUP_STRUCT *groups;
	int groupNum;
	LIST_STRUCT *lists;
	int listNum;
	char inputFileName[1000], outputFileName[1000];
	double maxPercentile;
	
	//Parse the command line
	if (argc == 1)
	{
		PrintCommandUsage(argv[0]);
		return -1;
	}
	
	inputFileName[0] = 0;
	outputFileName[0] = 0;
	maxPercentile = 0.25;
	
	for (i=2;i<argc;i++)
	{
		if (strcmp(argv[i-1], "-i")==0)
		{
			strcpy(inputFileName, argv[i]);
		}
		if (strcmp(argv[i-1], "-o")==0)
		{
			strcpy(outputFileName, argv[i]);
		}
		if (strcmp(argv[i-1], "-p")==0)
		{
			maxPercentile = atof(argv[i]);
		}
	}
	
	if ((inputFileName[0]==0)||(outputFileName[0]==0))
	{
		printf("Command error!\n");
		PrintCommandUsage(argv[0]);
		return -1;
	}
	
	if ((maxPercentile>1.0)||(maxPercentile<0.0))
	{
		printf("maxPercentile should be within 0.0 and 1.0\n");
		printf("program exit!\n");
		return -1;
	}
	
	groups = (GROUP_STRUCT *)malloc(MAX_GROUP_NUM*sizeof(GROUP_STRUCT));
	lists = (LIST_STRUCT *)malloc(MAX_LIST_NUM*sizeof(LIST_STRUCT));
	assert(groups!=NULL);
	assert(lists!=NULL);
	
	printf("reading input file...");
	
	flag = ReadFile(inputFileName, groups, MAX_GROUP_NUM, &groupNum, lists, MAX_LIST_NUM, &listNum);
	
	if (flag<=0)
	{
		printf("\nfailed.\n");
		printf("program exit!\n");
		
		return -1;
	}
	else
	{
		printf("done.\n");
	}
	
	printf("computing lo-values for each group...");
	
	if (ProcessGroups(groups, groupNum, lists, listNum, maxPercentile)<=0)
	{
		printf("\nfailed.\n");
		printf("program exit!\n");
		
		return -1;
	}
	else
	{
		printf("done.\n");
	}
	
	printf("computing false discovery rate...");
	
	if (ComputeFDR(groups, groupNum, maxPercentile, RAND_PASS_NUM*groupNum)<=0)
	{
		printf("\nfailed.\n");
		printf("program exit!\n");
		
		return -1;
	}
	else
	{
		printf("done.\n");
	}
	
	printf("save to output file...");
	
	if (SaveGroupInfo(outputFileName, groups, groupNum)<=0)
	{
		printf("\nfailed.\n");
		printf("program exit!\n");
		
		return -1;
	}
	else
	{
		printf("done.\n");
	}
	
	printf("finished.\n");
	
	free(groups);
	
	for (i=0;i<listNum;i++)
	{
		free(lists[i].values);
	}
	free(lists);

	return 0;

}

//print the usage of Command
void PrintCommandUsage(const char *command)
{
	//print the options of the command
	printf("%s - Robust Rank Aggreation.\n", command);
	printf("usage:\n");
	printf("-i <input data file>. Format: <item id> <group id> <list id> <value>\n");
	printf("-o <output file>. Format: <group id> <number of items in the group> <lo-value> <false discovery rate>\n");
	printf("-p <maximum percentile>. RRA only consider the items with percentile smaller than this parameter. Default=0.25\n");
	printf("example:\n");
	printf("%s -i input.txt -o output.txt -p 0.25 \n", command);
	
}

//Read input file. File Format: <item id> <group id> <list id> <value>. Return 1 if success, -1 if failure

int ReadFile(char *fileName, GROUP_STRUCT *groups, int maxGroupNum, int *groupNum, LIST_STRUCT *lists, int maxListNum, int *listNum)
{
	FILE *fh;
	int i,j;
	char **words, *tmpS;
	int wordNum;
	int totalItemNum;
	int tmpGroupNum, tmpListNum;
	char tmpGroupName[MAX_NAME_LEN], tmpListName[MAX_NAME_LEN], tmpItemName[MAX_NAME_LEN];
	double tmpValue;
	
	words = AllocWords(255, MAX_NAME_LEN+1);
	
	assert(words!=NULL);
	
	if (words == NULL)
	{
		return -1;
	}
	
	tmpS = (char *)malloc(255*(MAX_NAME_LEN+1)*sizeof(char));
	
	assert(tmpS!=NULL);
	
	fh = (FILE *)fopen(fileName, "r");
	
	if (!fh)
	{
		printf("Cannot open file %s\n", fileName);
		return -1;
	}
	
	//Read the header row to get the sample number
	fgets(tmpS, 255*(MAX_NAME_LEN+1)*sizeof(char), fh);
	
	wordNum = StringToWords(words, tmpS, MAX_NAME_LEN+1, 255, " \t\r\n\v\f");
	
	assert(wordNum == 4);
	
	if (wordNum != 4)
	{
		printf("Input file format: <item id> <group id> <list id> <value>\n");
		return -1;
	}
	
	//read records of items
	
	tmpGroupNum = 0;
	tmpListNum = 0;
	totalItemNum = 0;
	
	fgets(tmpS, 255*(MAX_NAME_LEN+1)*sizeof(char), fh);
	wordNum = StringToWords(words, tmpS, MAX_NAME_LEN+1, 255, " \t\r\n\v\f");
	
	while ((wordNum==4)&&(!feof(fh)))
	{
		strcpy(tmpItemName, words[0]);
		strcpy(tmpGroupName, words[1]);
		strcpy(tmpListName, words[2]);
		
		for (i=0;i<tmpGroupNum;i++)
		{
			if (!strcmp(tmpGroupName, groups[i].name))
			{
				break;
			} 
		}
		
		if (i>=tmpGroupNum)
		{
			strcpy(groups[tmpGroupNum].name, tmpGroupName);
			groups[tmpGroupNum].itemNum = 1;
			tmpGroupNum ++;
			if (tmpGroupNum >= maxGroupNum)
			{
				printf("too many groups. maxGroupNum = %d\n", maxGroupNum);
				return -1;
			}
		}
		else
		{
			groups[i].itemNum++;
		}
		
		for (i=0;i<tmpListNum;i++)
		{
			if (!strcmp(tmpListName, lists[i].name))
			{
				break;
			}
		}
		
		if (i>=tmpListNum)
		{
			strcpy(lists[tmpListNum].name, tmpListName);
			lists[tmpListNum].itemNum = 1;
			tmpListNum ++;
			if (tmpListNum >= maxListNum)
			{
				printf("too many lists. maxListNum = %d\n", maxListNum);
				return -1;
			}
		}
		else
		{
			lists[i].itemNum++;
		}
		
		totalItemNum++;
		
		fgets(tmpS, 255*(MAX_NAME_LEN+1)*sizeof(char), fh);
		wordNum = StringToWords(words, tmpS, MAX_NAME_LEN+1, 255, " \t\r\n\v\f");
	}
	
	fclose(fh);
	
	for (i=0;i<tmpGroupNum;i++)
	{
		groups[i].items = (ITEM_STRUCT *)malloc(groups[i].itemNum*sizeof(ITEM_STRUCT));
		groups[i].itemNum = 0;
	}
	
	for (i=0;i<tmpListNum;i++)
	{
		lists[i].values = (double *)malloc(lists[i].itemNum*sizeof(double));
		lists[i].itemNum = 0;
	}
	
	fh = (FILE *)fopen(fileName, "r");
	
	if (!fh)
	{
		printf("Cannot open file %s\n", fileName);
		return -1;
	}
	
	//Read the header row to get the sample number
	fgets(tmpS, 255*(MAX_NAME_LEN+1)*sizeof(char), fh);
	
	//read records of items
	
	fgets(tmpS, 255*(MAX_NAME_LEN+1)*sizeof(char), fh);
	wordNum = StringToWords(words, tmpS, MAX_NAME_LEN+1, 255, " \t\r\n\v\f");
	
	while ((wordNum==4)&&(!feof(fh)))
	{
		strcpy(tmpItemName, words[0]);
		strcpy(tmpGroupName, words[1]);
		strcpy(tmpListName, words[2]);
		tmpValue = atof(words[3]);
		
		for (i=0;i<tmpGroupNum;i++)
		{
			if (!strcmp(tmpGroupName, groups[i].name))
			{
				break;
			} 
		}
		
		for (j=0;j<tmpListNum;j++)
		{
			if (!strcmp(tmpListName, lists[j].name))
			{
				break;
			}
		}
		
		assert(i<tmpGroupNum);
		assert(j<tmpListNum);
		
		strcpy(groups[i].items[groups[i].itemNum].name,tmpItemName);
		groups[i].items[groups[i].itemNum].value = tmpValue;
		groups[i].items[groups[i].itemNum].listIndex = j;
		groups[i].itemNum ++;
		
		lists[j].values[lists[j].itemNum] = tmpValue;
		lists[j].itemNum ++;
		
		fgets(tmpS, 255*(MAX_NAME_LEN+1)*sizeof(char), fh);
		wordNum = StringToWords(words, tmpS, MAX_NAME_LEN+1, 255, " \t\r\n\v\f");
	}
	
	fclose(fh);
	
	printf("%d items\n%d groups\n%d lists\n", totalItemNum, tmpGroupNum, tmpListNum);
	
	*groupNum = tmpGroupNum;
	*listNum = tmpListNum;
	
	FreeWords(words, 255);
	free(tmpS);
	
	return totalItemNum;
	
}

//Save group information to output file. Format <group id> <number of items in the group> <lo-value> <false discovery rate>
int SaveGroupInfo(char *fileName, GROUP_STRUCT *groups, int groupNum)
{
	FILE *fh;
	int i;
	
	fh = (FILE *)fopen(fileName, "w");
	
	if (!fh)
	{
		printf("Cannot open %s.\n", fileName);
		return -1;
	}
	
	fprintf(fh, "group_id\t#_items_in_group\tlo_value\tFDR\n");
	
	for (i=0;i<groupNum;i++)
	{
		fprintf(fh, "%s\t%d\t%10.4e\t%f\n", groups[i].name, groups[i].itemNum, groups[i].loValue, groups[i].fdr);
	}
	
	fclose(fh);
	
	return 1;
}

//Process groups by computing percentiles for each item and lo-values for each group
int ProcessGroups(GROUP_STRUCT *groups, int groupNum, LIST_STRUCT *lists, int listNum, double maxPercentile)
{
	int i,j;
	int listIndex, index1, index2;
	int maxItemPerGroup;
	double *tmpF;
	
	maxItemPerGroup = 0;
	
	for (i=0;i<groupNum;i++)
	{
		if (groups[i].itemNum>maxItemPerGroup)
		{
			maxItemPerGroup = groups[i].itemNum;
		}
	}
	
	assert(maxItemPerGroup>0);
	
	tmpF = (double *)malloc(maxItemPerGroup*sizeof(double));
	
	for (i=0;i<listNum;i++)
	{
		QuicksortF(lists[i].values, 0, lists[i].itemNum-1);
	}
	
	for (i=0;i<groupNum;i++)
	{
		//Compute percentile for each item
		
		for (j=0;j<groups[i].itemNum;j++)
		{
			listIndex = groups[i].items[j].listIndex;
			
			index1 = bTreeSearchingF(groups[i].items[j].value-0.000000001, lists[listIndex].values, 0, lists[listIndex].itemNum-1);
			index2 = bTreeSearchingF(groups[i].items[j].value+0.000000001, lists[listIndex].values, 0, lists[listIndex].itemNum-1);
			
			groups[i].items[j].percentile = ((double)index1+index2+1)/(lists[listIndex].itemNum*2);
			tmpF[j] = groups[i].items[j].percentile;
		}
		
		ComputeLoValue(tmpF, groups[i].itemNum, &(groups[i].loValue), maxPercentile);
	}
	
	return 1;
}

//Compute lo-value based on an array of percentiles. Return 1 if success, 0 if failure
int ComputeLoValue(double *percentiles,     //array of percentiles
				   int num,                 //length of array
				   double *loValue,         //pointer to the output lo-value
				   double maxPercentile)   //maximum percentile, computation stops when maximum percentile is reached
{
	int i;
	double *tmpArray;
	double tmpLoValue, tmpF;
	
	assert(num>0);
	
	tmpArray = (double *)malloc(num*sizeof(double));
	
	if (!tmpArray)
	{
		return -1;
	}
	
	memcpy(tmpArray, percentiles, num*sizeof(double));
	
	QuicksortF(tmpArray, 0, num-1);
	
	tmpLoValue = 1.0;
	
	for (i=0;i<num;i++)
	{
		if ((tmpArray[i]>maxPercentile)&&(i>0))
		{
			break;
		}
		tmpF = BetaNoncentralCdf((double)(i+1),(double)(num-i),0.0,tmpArray[i],CDF_MAX_ERROR);
		if (tmpF<tmpLoValue)
		{
			tmpLoValue = tmpF;
		}
	}
	
	*loValue = tmpLoValue;

	free(tmpArray);
	
	return 0;
	
}

//Compute False Discovery Rate based on uniform distribution
int ComputeFDR(GROUP_STRUCT *groups, int groupNum, double maxPercentile, int numOfRandPass)
{
	int i,j,k;
	double *tmpPercentile;
	int maxItemNum = 0;
	int scanPass = numOfRandPass/groupNum+1;
	double *randLoValue;
	int randLoValueNum;
	
	for (i=0;i<groupNum;i++)
	{
		if (groups[i].itemNum>maxItemNum)
		{
			maxItemNum = groups[i].itemNum;
		}
	}
	
	assert(maxItemNum>0);
	
	tmpPercentile = (double *)malloc(maxItemNum*sizeof(double));
	
	randLoValueNum = groupNum*scanPass;
	
	assert(randLoValueNum>0);
	
	randLoValue = (double *)malloc(randLoValueNum*sizeof(double));
	
	randLoValueNum = 0;
	
	PlantSeeds(123456);
	
	for (i=0;i<scanPass;i++)
	{
		for (j=0;j<groupNum;j++)
		{
			for (k=0;k<groups[j].itemNum;k++)
			{
				tmpPercentile[k] = Uniform(0.0, 1.0);
			}
			
			ComputeLoValue(tmpPercentile, groups[j].itemNum,&(randLoValue[randLoValueNum]), maxPercentile);
			
			randLoValueNum++;
		}
	}
	
	QuicksortF(randLoValue, 0, randLoValueNum-1);
						  
	QuickSortGroupByLoValue(groups, 0, groupNum-1);
	
	for (i=0;i<groupNum;i++)
	{
		groups[i].fdr = (double)(bTreeSearchingF(groups[i].loValue-0.000000001, randLoValue, 0, randLoValueNum-1)
								 +bTreeSearchingF(groups[i].loValue+0.000000001, randLoValue, 0, randLoValueNum-1)+1)
								/2/randLoValueNum/((double)i+0.5)*groupNum;
	}
	
	if (groups[groupNum-1].fdr>1.0)
	{
		groups[groupNum-1].fdr = 1.0;
	}
	
	for (i=groupNum-2;i>=0;i--)
	{
		if (groups[i].fdr>groups[i+1].fdr)
		{
			groups[i].fdr = groups[i+1].fdr;
		}
	}
	
	free(tmpPercentile);
	free(randLoValue);
	
	return 1;
}

//QuickSort groups by loValue
void QuickSortGroupByLoValue(GROUP_STRUCT *groups, int lo, int hi)
{
	int i=lo, j=hi;
	GROUP_STRUCT tmpGroup;
	double x=groups[(lo+hi)/2].loValue;
	
	if (hi<lo)
	{
		return;
	}
	
    //  partition
    while (i<=j)
    {    
		while ((groups[i].loValue<x)&&(i<=j))
		{
			i++;
		}
		while ((groups[j].loValue>x)&&(i<=j))
		{
			j--;
		}
        if (i<=j)
        {
			memcpy(&tmpGroup,groups+i,sizeof(GROUP_STRUCT));
			memcpy(groups+i,groups+j,sizeof(GROUP_STRUCT));
			memcpy(groups+j,&tmpGroup,sizeof(GROUP_STRUCT));
            i++; j--;
        }
    } 
	
    //  recursion
    if (lo<j) QuickSortGroupByLoValue(groups, lo, j);
    if (i<hi) QuickSortGroupByLoValue(groups, i, hi);
	
}
