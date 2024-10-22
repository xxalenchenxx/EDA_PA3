/* -----------  PA3  ----------------
    Chen-Yu-Chia
    2024-10-22
------------------------------------- */
/* --------------------------------------------------------------------------
   Header file used in bookshelf_IO.c 
----------------------------------------------------------------------------*/

#ifndef _BOOKSHELF_IO_H_
#define _BOOKSHELF_IO_H_

#define BUFFERSIZE 300
#define STRINGLEN 40

    /* -----------------------------------------------------------------------------
        Reads the .nodes file and creates a hash linking cell name to cell index for 
        all the nodes in the circuit (movable nodes + fixed nodes + I/O pads)

        creates extern vars:
            cellName[i]     (i = 1..movableNodes + numTerminals)
    -------------------------------------------------------------------------------- */
    extern void createHash(char benchmarkPath[], char nodesFile[]);
    extern void freeHash();

    /* -----------------------------------------------------------------------------
        Reads the .aux file to get the other file names
  
        creates extern vars:
            nodesFile[], netsFile[], wtsFile[], sclFile[], plFile[], benchmarkName[]
    -------------------------------------------------------------------------------- */
    extern void readAuxFile(char benchmarkPath[], char auxFile[]);

    /* -----------------------------------------------------------------------------
        Reads the .nodes file to get cell widths and heights
        
        creates extern vars: 
            movableNodes, numTerminals, averageCellWidth, cellWidth[], cellHeight[]
    -------------------------------------------------------------------------------- */
    extern void readNodesFile(char benchmarkPath[], char nodesFile[]);

    /* -----------------------------------------------------------------------------
        Reads the .nets file to get the netlist information
   
        creates extern vars: 
            numNets, numPins, netlist[], netlistIndex[], xPinOffset[], yPinOffset[]
    -------------------------------------------------------------------------------- */
    extern void readNetsFile(char benchmarkPath[], char netsFile[]);

    /* -----------------------------------------------------------------------------
        Reads the .pl file to get coordinates of all the nodes and the placement 
        boundary based on the position of the I/O pads
  
        creates extern vars:
            xCellCoord[], yCellCoord[]
    -------------------------------------------------------------------------------- */
    extern void readPlFile(char benchmarkPath[], char plFile[]);
    
    /* -----------------------------------------------------------------------------
        Reads the .scl file to get placement (core) region information
  
        creates extern vars:
            siteOriginX, siteEndX, siteOriginY, siteEndY, siteWidth, siteSpacing, 
            numRows, coreRowHeight, coreWidth, coreHeight 
    -------------------------------------------------------------------------------- */
    extern void readSclFile(char benchmarkPath[], char sclFile[]);
    
    /* -----------------------------------------------------------------------------
        writes out a bookshelf format .pl file
    -------------------------------------------------------------------------------- */
    extern void writePlFile(char outputDir[], char benchmarkName[], float xCoord[], float yCoord[]);

    /* -----------------------------------------------------------------------------
        print out a bookshelf format .pl file to show in python
    -------------------------------------------------------------------------------- */
    extern void write_python_File(float *xCellCoord,float *yCellCoord);


    /* -----------------------------------------------------------------------------
        output results into output directory
    -------------------------------------------------------------------------------- */
    extern void clear_output_dir(char outputDir[]);
    extern void write_aux_File(char outputDir[], char benchmarkName[]);

    //write_other_File: ".nodes", ".nets", ".wts",  ".scl"
    extern void write_other_File(char inputDir[], char outputDir[], char benchmarkName[]);
    /*--------------  Extern Variables  ------------------*/

    extern char **cellName;
    
    extern char nodesFile[BUFFERSIZE], netsFile[BUFFERSIZE], wtsFile[BUFFERSIZE];
    extern char sclFile[BUFFERSIZE], plFile[BUFFERSIZE], benchmarkName[BUFFERSIZE];

    extern int movableNodes, numTerminals;
    extern float averageCellWidth, *cellWidth, *cellHeight; 
    
    extern int numNets, numPins, *netlist, *netlistIndex;
    extern float *xPinOffset, *yPinOffset;
       
    extern float *xCellCoord, *yCellCoord, minX, maxX, minY, maxY;
    extern char **cell_fixedType; 
    extern int *areaArrayIO, numAreaArrayIO;

    extern int numRows, numRowBlockages;
    extern float siteOriginY, siteEndY, coreHeight;
    extern float siteOriginX, siteEndX, coreWidth;
    extern float siteWidth, siteSpacing, coreRowHeight;
    extern float *rowOriginX, *rowEndX;
    extern float *xRowBlockage, *yRowBlockage, *widthRowBlockage,*siteSpacingRow;
    extern int numNodes;

#endif /* _BOOKSHELF_IO_H_*/ 
