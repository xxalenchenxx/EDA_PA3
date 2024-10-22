#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <cmath>
#include <ctime>
#include <vector> 
#include <limits>
#include <cstdlib>
#include <iostream>
#include <algorithm>        
#include "bookshelf_IO.h"
#include "memAlloc.h"

struct Rect {
    float xStart, xEnd;
    float yStart, yEnd;
    int index; 
};


float clamp(float val, float minVal, float maxVal) {
    return std::max(minVal, std::min(maxVal, val));
}


bool compareRect(const Rect &a, const Rect &b) {
    if (a.yStart != b.yStart)
        return a.yStart < b.yStart; // sort y by up-down
    else
        return a.xStart < b.xStart; // if y is same, I sort X by left-right
}

void check_overlap(std::vector<float> xCellCoord_p,std::vector<float> yCellCoord_p){
    // 检查是否有重叠的矩形
    bool hasOverlap = false;
    for (int i = 1; i <= numNodes; ++i) {
        for (int j = i + 1; j <= numNodes; ++j) {
            // 获取矩形 i 的范围
            float xStart_i = xCellCoord_p[i];
            float xEnd_i = xStart_i + cellWidth[i];
            float yStart_i = yCellCoord_p[i];
            float yEnd_i = yStart_i + cellHeight[i];

            // 获取矩形 j 的范围
            float xStart_j = xCellCoord_p[j];
            float xEnd_j = xStart_j + cellWidth[j];
            float yStart_j = yCellCoord_p[j];
            float yEnd_j = yStart_j + cellHeight[j];

            // 检查 X 轴上的重叠，边界相等不算重叠
            bool xOverlap = xStart_i < xEnd_j && xEnd_i > xStart_j;
            if (xEnd_i == xStart_j || xEnd_j == xStart_i) {
                xOverlap = false;
            }

            // 检查 Y 轴上的重叠，边界相等不算重叠
            bool yOverlap = yStart_i < yEnd_j && yEnd_i > yStart_j;
            if (yEnd_i == yStart_j || yEnd_j == yStart_i) {
                yOverlap = false;
            }

            if (xOverlap && yOverlap) {
                hasOverlap = true;
                // std::cout << "Overlap detected between Rect " << cellName[i] <<" ("<<xStart_i<<" , "<<yStart_i<< " ) and Rect " 
                //           << cellName[j] <<" ("<<xStart_j<<" , "<<yStart_j<< " )" << std::endl;
            }
        }
    }

    if (!hasOverlap) {
        std::cout << "[CORRECT]No overlaps" << std::endl;
    }
    return;
}


void check_on_site(std::vector<float> xCellCoord_p,std::vector<float> yCellCoord_p,std::vector<float> siteSpacing,std::vector<float> subrowOrigin){
//檢查是否符合 site的座標
    bool has_non_site = false;
    for (int i = 1; i <= numNodes; ++i) {
        float xStart_i = xCellCoord_p[i];
        float yStart_i = yCellCoord_p[i];

        // 確定 cell 所在的 row
        int row = static_cast<int>((yStart_i - siteOriginY) / coreRowHeight);
        // 檢查 row 索引是否合法
        if (row < 0 || row >= numRows) {
            std::cout << "Error: Cell " << i << " is in invalid row." << std::endl;
            has_non_site = true;
            break;
        }

        // 獲取該 row 的 Sitespacing 和 SubrowOrigin
        float spacing = siteSpacing[row];
        float origin = subrowOrigin[row];

        // 計算 x 座標相對於 SubrowOrigin 的偏移量
        float offset = xStart_i - origin;

        // 檢查偏移量是否為 spacing 的整數倍
        float quotient = offset / spacing;
        int quotient_int = static_cast<int>(std::round(quotient));
        float diff = std::fabs(offset - quotient_int * spacing);
        const float epsilon = 1e-6; // 允許的微小誤差

        if (diff > epsilon) {
            std::cout << "Cell " << i << " xStart_i: " << xStart_i << " is not aligned to site grid in row " << row << std::endl;
            has_non_site = true;
            break;
        }
    }
    if (!has_non_site) {
        std::cout << "[CORRECT] All positions are on site grid" << std::endl;
    }

}


bool checkOverlap(const std::vector<float> &xCellCoord_p, const std::vector<float> &yCellCoord_p, int numNodes) {
    for (int i = 1; i <= numNodes; ++i) {
        for (int j = i + 1; j <= numNodes; ++j) {
            // 获取矩形 i 的范围
            float xStart_i = xCellCoord_p[i];
            float xEnd_i = xStart_i + cellWidth[i];
            float yStart_i = yCellCoord_p[i];
            float yEnd_i = yStart_i + cellHeight[i];

            // 获取矩形 j 的范围
            float xStart_j = xCellCoord_p[j];
            float xEnd_j = xStart_j + cellWidth[j];
            float yStart_j = yCellCoord_p[j];
            float yEnd_j = yStart_j + cellHeight[j];

            // 检查 X 轴上的重叠，边界相等不算重叠
            bool xOverlap = xStart_i < xEnd_j && xEnd_i > xStart_j;
            if (xEnd_i == xStart_j || xEnd_j == xStart_i) {
                xOverlap = false;
            }

            // 检查 Y 轴上的重叠，边界相等不算重叠
            bool yOverlap = yStart_i < yEnd_j && yEnd_i > yStart_j;
            if (yEnd_i == yStart_j || yEnd_j == yStart_i) {
                yOverlap = false;
            }

            if (xOverlap && yOverlap) {
                return true; // 有重叠
            }
        }
    }
    return false; // 无重叠
}


void calculateDisplacement(const std::vector<float> &xCellCoord_p, const std::vector<float> &yCellCoord_p, int numNodes, float &total_cost) {
    total_cost = 0.0f;
    for (int i = 1; i <= numNodes; ++i) {
        float ab_x = xCellCoord_p[i] - xCellCoord[i];
        float ab_y = yCellCoord_p[i] - yCellCoord[i];
        float cost = std::fabs(ab_x) + std::fabs(ab_y);
        total_cost += cost;

    }
}



void simulatedAnnealing(std::vector<float> &xCellCoord_p, std::vector<float> &yCellCoord_p,
                        float *cellWidth, float *cellHeight,
                        const std::vector<float> &siteSpacing, const std::vector<float> &subrowOrigin,
                        int numNodes, int numRows, float siteOriginY, float coreRowHeight) {
    // 参数设置
    float T = 1024.0f; // intial T
    float T_min = 1e-1f; // minimal T
    float alpha = 0.997f; // coe
    int maxIterations = 2048; // max iteration

    // 构建宽度到cell索引的映射
    std::map<float, std::vector<int>> widthToCells;
    for (int i = 1; i <= numNodes; ++i) {
        widthToCells[cellWidth[i]].push_back(i);
    }

    std::vector<float> currentX = xCellCoord_p;
    std::vector<float> currentY = yCellCoord_p;


    float currentTotalCost, currentMaxCost;
    calculateDisplacement(currentX, currentY, numNodes, currentTotalCost);
    float historyTotalCost=1e10;
    std::srand(static_cast<unsigned int>(std::time(nullptr))); // random


    //start SA
    while (T > T_min) {
        // std::cout<<"T:"<<T<<std::endl;
        for (int iter = 0; iter < maxIterations; ++iter) {
            
            auto it = widthToCells.begin();
            std::advance(it, std::rand() % widthToCells.size());
            float width = it->first;
            auto &cells = it->second;

            if (cells.size() < 2) continue;

            
            int idx1 = cells[std::rand() % cells.size()];
            int idx2 = cells[std::rand() % cells.size()];
            float dist=std::fabs(currentX[idx1]-currentX[idx2])+std::fabs(currentY[idx1]-currentY[idx2]);
            if (idx1 == idx2 || cellWidth[idx1]!=cellWidth[idx2]) continue;
            // else if(dist>10000){
            //     continue;
            // }

            
            std::swap(currentX[idx1], currentX[idx2]);
            std::swap(currentY[idx1], currentY[idx2]);

            

            // new cost
            float newTotalCost;
            calculateDisplacement(currentX, currentY, numNodes, newTotalCost);

            if (newTotalCost < currentTotalCost) {
                currentTotalCost = newTotalCost;
            } else {
                float delta = newTotalCost - currentTotalCost;
                float probability = std::exp(-delta / T);
                if ((std::rand() / (float)RAND_MAX) < probability) {

                    currentTotalCost = newTotalCost;
                    // printf("[BAD]Total displacement: %.1f\n",currentTotalCost);
                } else {

                    std::swap(currentX[idx1], currentX[idx2]);
                    std::swap(currentY[idx1], currentY[idx2]);
                }
            }
            if(currentTotalCost<historyTotalCost){
                // printf("[GOOD]Total displacement: %.1f\n",currentTotalCost);
                historyTotalCost=currentTotalCost;
                xCellCoord_p = currentX;
                yCellCoord_p = currentY;
            }

        }

        // if(T <100.0){
        //     T *= (alpha+0.20);
        // }else{
            T *= alpha;
        // }
        
    }

}

int main (int argc, char *argv[])
{
    char auxFile[BUFFERSIZE], benchmarkPath[BUFFERSIZE], outputPath[BUFFERSIZE];

    if(argc != 4) {
        printf("Usage: %s <benchmark_dir> <aux_file> <placement_file>\n",
               argv[0]);
        printf("    <benchmark_dir> is the benchmark file directory.\n");
        printf("    <aux_file> is the bookshelf format auxiliary file");
        printf("  (assume in <benchmark_dir>).\n");
        printf("    <placement_file> is the placement file");
        printf(" (assume in current directory).\n");
        exit(1);
    }    

    strcpy(benchmarkPath, argv[1]);
    strcpy(auxFile, argv[2]);
    strcpy(outputPath, argv[3]);

    readAuxFile(benchmarkPath, auxFile);
    createHash(benchmarkPath, nodesFile);
    readNodesFile(benchmarkPath, nodesFile);

    readNetsFile(benchmarkPath, netsFile);
    readPlFile(benchmarkPath, plFile);
    readSclFile(benchmarkPath, sclFile);

    /* -----------------------------------------------------------------------------
        Modified Tetris Legalization Algorithm
    -------------------------------------------------------------------------------- */

    // define placed position (left down)
    std::vector<float> xCellCoord_p(numNodes + 1);
    std::vector<float> yCellCoord_p(numNodes + 1);

    // construct rect
    std::vector<Rect> rects;
    for (int i = 1; i <= numNodes; ++i) {
        Rect rect;
        rect.xStart = xCellCoord[i];
        rect.yStart = yCellCoord[i];
        rect.xEnd = rect.xStart + cellWidth[i];
        rect.yEnd = rect.yStart + cellHeight[i];
        rect.index = i;
        rects.push_back(rect);
    }


    std::sort(rects.begin(), rects.end(), compareRect);


    std::vector<std::vector<std::pair<float, float>>> availableIntervals(numRows);
    // store row of Sitespacing & SubrowOrigin
    std::vector<float> siteSpacing(numRows);
    std::vector<float> subrowOrigin(numRows);
    std::vector<int> numSites(numRows);


    for (int row = 0; row < numRows; ++row) {
        availableIntervals[row].push_back({rowOriginX[row], rowEndX[row]});
        siteSpacing[row] = siteSpacingRow[row+1]; 
        subrowOrigin[row] = rowOriginX[row]; // SubrowOrigin 是 row 的起始 x 座標
        numSites[row] = (int)((rowEndX[row] - rowOriginX[row]) / siteSpacing[row]);
    }

    // Modified Tetris Legalization
    for (auto &rect : rects) {
        float width = cellWidth[rect.index];
        float height = cellHeight[rect.index];

        // 存储最佳位置信息
        bool placed = false;
        float minDistance = std::numeric_limits<float>::max();
        int bestRow = -1;
        float bestXStart = rect.xStart;
        float bestYStart = rect.yStart;

        // 遍历所有 row，寻找最适合的位置
        for (int row = 0; row < numRows; ++row) {
            float yStart = siteOriginY + row * coreRowHeight;

            // 确保矩形高度不超过 row 的高度
            if (height > coreRowHeight) {
                continue; // 无法放置在该 row
            }

            // 获取该 row 的 Sitespacing 和 SubrowOrigin
            float spacing = siteSpacing[row]; //66
            float origin = subrowOrigin[row]; 
            int sites = numSites[row]; //66

            // 在该 row 中遍历可用区间
            auto &intervals = availableIntervals[row];
            for (auto it = intervals.begin(); it != intervals.end(); ++it) {
                float intervalStart = it->first;
                float intervalEnd = it->second;

                // 计算区间内的起始和结束 site 索引
                int startSiteIndex = std::max(0, (int)std::ceil((intervalStart - origin) / spacing));
                int endSiteIndex = std::min(sites, (int)((intervalEnd - origin) / spacing));

                // 遍历该区间内的可能 site 位置
                for (int siteIndex = startSiteIndex; siteIndex < endSiteIndex; ++siteIndex) {
                    float xSite = origin + siteIndex * spacing;

                    // 确保矩形能放入区间
                    if (xSite + width >= intervalEnd ) { // 允许小的误差
                        break; // 不能再放置更多 site
                    }

                    // 计算与原始位置的距离
                    float distance = std::fabs(xSite - rect.xStart) + std::fabs(yStart - rect.yStart);

                    if (distance < minDistance) {
                        minDistance = distance;
                        bestXStart = xSite;
                        bestYStart = yStart;
                        bestRow = row;
                        placed = true;
                    }
                }
            }
        }

        if (placed) {
            // 更新最佳 row 的可用区间
            float xPlacedStart = bestXStart;
            float xPlacedEnd = xPlacedStart + width;

            auto &intervals = availableIntervals[bestRow];
            std::vector<std::pair<float, float>> newIntervals;

            for (auto &interval : intervals) {
                if (interval.second <= xPlacedStart || interval.first >= xPlacedEnd) {
                    // 不重叠，保留
                    newIntervals.push_back(interval);
                } else {
                    // 可能需要拆分区间
                    if (interval.first < xPlacedStart) {
                        newIntervals.push_back({interval.first, xPlacedStart});
                    }
                    if (interval.second > xPlacedEnd) {
                        newIntervals.push_back({xPlacedEnd, interval.second});
                    }
                }
            }
            intervals = newIntervals;

            // 记录新位置
            xCellCoord_p[rect.index] = bestXStart;
            yCellCoord_p[rect.index] = bestYStart;
        } else {
            // 无法放置，保持原始位置并输出警告
            xCellCoord_p[rect.index] = rect.xStart;
            yCellCoord_p[rect.index] = rect.yStart;
            std::cout << "Warning: Rect " << cellName[rect.index] << " could not be placed." <<
                " ("<<xCellCoord_p[rect.index]<<" , "<<yCellCoord_p[rect.index]<< " )"<< std::endl;
        }
    }


    float total_cost=0.0;
    float max_cost=0.0;
    for(int i=1;i<xCellCoord_p.size();i++){
        float ab_x=xCellCoord_p[i]-xCellCoord[i];
        float ab_y=yCellCoord_p[i]-yCellCoord[i];
        float cost= std::fabs(ab_x)+std::fabs(ab_y);
        total_cost+=cost;
        if(max_cost<cost){
            max_cost=cost;
        }
    }
    printf("[before]Total displacement: %.1f\n",total_cost);
    printf("[before]Maximum displacement: %.1f\n",max_cost);

    /* -----------------------------------------------------------------------------
        simulated Annealing Algorithm 
    -------------------------------------------------------------------------------- */
    
    simulatedAnnealing(xCellCoord_p, yCellCoord_p, cellWidth, cellHeight, siteSpacing, subrowOrigin,
                       numNodes, numRows, siteOriginY, coreRowHeight);




    // visualize result
    write_python_File( xCellCoord_p.data(), yCellCoord_p.data());


    //Total cost & MAX cost
    total_cost=0.0;
    max_cost=0.0;
    for(int i=1;i<xCellCoord_p.size();i++){
        float ab_x=xCellCoord_p[i]-xCellCoord[i];
        float ab_y=yCellCoord_p[i]-yCellCoord[i];
        float cost= std::fabs(ab_x)+std::fabs(ab_y);
        total_cost+=cost;
        if(max_cost<cost){
            max_cost=cost;
        }
    }
    printf("[after]Total displacement: %.1f\n",total_cost);
    printf("[after]Maximum displacement: %.1f\n",max_cost);

    // check_overlap
    check_overlap(xCellCoord_p,yCellCoord_p);
    //檢查是否符合 site的座標
    check_on_site(xCellCoord_p,yCellCoord_p,siteSpacing,subrowOrigin);

    freeHash();
    return 0;
}
