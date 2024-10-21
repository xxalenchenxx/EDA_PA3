#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "bookshelf_IO.h"
#include "memAlloc.h"
#include <vector>
#include <cmath>
#include <math.h>
#include <iostream>
#include <limits>
#include <cstdlib> // 用於 rand()
#include <ctime>   // 用於時間種子

// 定義一個簡單的結構來表示 rect 的範圍
struct Rect {
    float xStart, xEnd;
    float yStart, yEnd;
};

// 檢查兩個矩形是否重疊
bool checkOverlap(const Rect &a, const Rect &b) {
    // 檢查 X 軸和 Y 軸上的重疊
    bool xOverlap = !(a.xEnd <= b.xStart || b.xEnd <= a.xStart);
    bool yOverlap = !(a.yEnd <= b.yStart || b.yEnd <= a.yStart);
    return xOverlap && yOverlap;
}

// 將值限制在最小值和最大值之間
float clamp(float val, float minVal, float maxVal) {
    return std::max(minVal, std::min(maxVal, val));
}

int main (int argc, char *argv[])
{
    char auxFile[BUFFERSIZE], benchmarkPath[BUFFERSIZE], outputfile[BUFFERSIZE];

    if(argc != 4) {
        printf("Usage: %s <benchmark_dir> <aux_file> <placement_file>\n",
               argv[0]);
        printf("    <benchmark_dir> is the benchmark file directory.\n");
        printf("    <aux_file> is the bookshelf format auxiliary file");
        printf(" (assume in <benchmark_dir>).\n");
        printf("    <placement_file> is the placement file");
        printf(" (assume in current directory).\n");
        exit(1);
    }    

    strcpy(benchmarkPath, argv[1]);
    strcpy(auxFile, argv[2]);
    strcpy(outputfile, argv[3]);

    readAuxFile(benchmarkPath, auxFile);
    createHash(benchmarkPath, nodesFile);
    readNodesFile(benchmarkPath, nodesFile);

    readNetsFile(benchmarkPath, netsFile);
    readPlFile(benchmarkPath, plFile);
    readSclFile(benchmarkPath, sclFile);

    // 初始化隨機數種子
    std::srand(static_cast<unsigned int>(std::time(nullptr)));

    /* -----------------------------------------------------------------------------
        Placement Algorithm
    -------------------------------------------------------------------------------- */

    float* xCellCoord_p = vector(1,numNodes);
    float* yCellCoord_p = vector(1,numNodes);

    // 創建已放置的矩形列表
    std::vector<Rect> placedRects;

    for (int i = 1; i <= numNodes; ++i) {
        float originalX = xCellCoord[i];
        float originalY = yCellCoord[i];
        float width = cellWidth[i];
        float height = cellHeight[i];

        // 初始搜索偏移量
        // float offset = coreRowHeight; // 可以根據需求調整初始偏移量
        bool placed = false;

        // 最大嘗試次數，防止無限循環
        int dist=0;
        while (!placed && dist<1000) {
            // 隨機選擇方向：0=左，1=右，2=上，3=下
            
            int direction = rand() % 8;
            float newXStart = originalX;
            float newYStart = std::ceil(originalY / cellHeight[i]) * cellHeight[i];
            // std::cout << "originalY: " << originalY << " ncellHeight[i] " << cellHeight[i] <<"newYStart "<<newYStart<< std::endl;
            float offset = (float)((int)(rand()%100)*0.001);
            if(dist){
                switch (direction) {
                    case 0: // 左
                        newXStart -= newXStart*offset;
                        break;
                    case 1: // 右
                        newXStart += newXStart*offset;
                        break;
                    case 2: // 上
                        newYStart += cellHeight[i]*dist;
                        break;
                    case 3: // 下
                        newYStart -= cellHeight[i]*dist;
                        break;

                    case 4: // 左上
                        newXStart -= newXStart*offset;
                        newYStart += cellHeight[i]*dist;
                        break;
                    case 5: // 右上
                        newXStart += newXStart*offset;
                        newYStart += cellHeight[i]*dist;
                        break;
                    case 6: // 右下
                        newXStart += newXStart*offset;
                        newYStart -= cellHeight[i]*dist;
                        break;
                    case 7: // 左下
                        newXStart -= newXStart*offset;
                        newYStart -= cellHeight[i]*dist;
                        break;
                    
                }
            }

            
            dist++;
            // 將矩形限制在範圍內
            newXStart = clamp(newXStart, siteOriginX, coreWidth - width);
            newYStart = clamp(newYStart, siteOriginY, siteEndY - height);

            // 計算矩形的終點
            float newXEnd = newXStart + width;
            float newYEnd = newYStart + height;

            Rect newRect = {newXStart, newXEnd, newYStart, newYEnd};
            // 檢查是否與已放置的矩形重疊
            bool overlaps = false;
            for (const auto &rect : placedRects) {
                if (checkOverlap(newRect, rect)) {
                    overlaps = true;
                    break;
                }
            }

            if (!overlaps) {
                // 無重疊，放置矩形
                xCellCoord_p[i] = newXStart;
                yCellCoord_p[i] = newYStart;
                placedRects.push_back(newRect);
                placed = true;
            }
            else{
                newXStart = originalX;
                newYStart = std::ceil(originalY / cellHeight[i]) * cellHeight[i];
            }
            
        }

        if (!placed) {
            // 無法放置矩形，輸出警告並保持原始位置
            xCellCoord_p[i] = originalX;
            yCellCoord_p[i] = originalY;
            // std::cout << "Warning: Rect " << i << " could not be placed without overlap after " << attempts << " attempts." << std::endl;
        }
    }

    // 輸出結果
    double cost=0.0;
    double cell_cost_max=0.0;
    for (int i = 1; i <= numNodes; ++i) {
        float x_cost=xCellCoord_p[i]-xCellCoord[i];
        float y_cost=yCellCoord_p[i]-yCellCoord[i];
        float cell_cost=std::fabs(x_cost)+std::fabs(y_cost);
        cost+=cell_cost;
        if(cell_cost>cell_cost_max){
            cell_cost_max=cell_cost;
        }
    }
    std::cout<<"placement cost: "<<cost<<std::endl;
    std::cout<<"cell max  cost: "<<cell_cost_max<<std::endl;
    write_python_File(xCellCoord_p, yCellCoord_p);

    freeHash();
    return 0;
}
