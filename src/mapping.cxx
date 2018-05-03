
#include "TdcMapping.hh"

#include <stdio.h>
int main()
{
  // printf(" TDC2PR %d PR2LEMO %d \n",TDC2PR[7],PR2LEMO[TDC2PR[7]]);
  // printf(" TDC2PR %d PR2LEMO %d \n",LEMO2PR[18],PR2TDC[LEMO2PR[18]]);

  // for (int itdc=0;itdc<32;itdc++)
  //   {
  //     int ilemo=-1;
  //     for (int j=0;j<32;j++)
  // 	if (LEMO2TDC[j]==itdc) ilemo=j;
  //     printf("%d,",ilemo);
 
  //   }
  // printf("\n ");



  //   for (int ipr=0;ipr<32;ipr++)
  //   {
  //     int ilemo=-1;
  //     if (PR2TDC[ipr]!=-1)
  // 	ilemo = TDC2LEMO[PR2TDC[ipr]];
  //     printf("%d,",ilemo);
 
  //   }
  // printf("\n");

  //  for (int itdc=0;itdc<32;itdc++)
  //   {
  //     int ilemo=-1;
  //     for (int j=0;j<32;j++)
  // 	if (PR2LEMO[j]==itdc) ilemo=j;
  //     printf("%d,",ilemo);
 
  //   }
  // printf("\n");



  //  for (int itdc=0;itdc<24;itdc++)
  //   {
  //     printf("TDC %d LEMO %d STRIP %d \n",itdc,TDC2LEMO[itdc],LEMO2STRIP[TDC2LEMO[itdc]]);
 
  //   }
  // printf("\n");
  // for (int istrip=71;istrip<86;istrip++)
  //   {
  //     printf("Strip %d \n ",istrip);
  //     for (int ilemo=0;ilemo<32;ilemo++)
  // 	if (LEMO2STRIP[ilemo]==istrip)
  // 	  {
  // 	    int pr = -1;
  // 	    if (LEMO2TDC[ilemo]>0)
  // 	      pr=TDC2PR[LEMO2TDC[ilemo]];
  // 	    printf("\t LEMO %d TDC %d  PR %d\n",ilemo,LEMO2TDC[ilemo],pr);
  // 	  }
  //   }
  // printf("\n");


   for (int ipr=0;ipr<32;ipr++)
    {
      int ilemo=-1;
      for (int j=0;j<32;j++)
	if (PR2TDC[j]==ipr && PR2TDC[j]>-1) ilemo=j;
      printf("TDC %d %d,\n",ipr,ilemo);
 
    }
  printf("\n");


}
