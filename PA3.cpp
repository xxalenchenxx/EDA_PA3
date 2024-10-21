#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "bookshelf_IO.h"
#include "memAlloc.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <limits>
#include <algorithm> // 用於 std::max 和 std::min

// 定義一個簡單的結構來表示 rect 的範圍
struct Rect {
    float xStart, xEnd;
    float yStart, yEnd;
    int index; // 紀錄矩形的索引
};

// 將值限制在最小值和最大值之間
float clamp(float val, float minVal, float maxVal) {
    return std::max(minVal, std::min(maxVal, val));
}

// 比較函數，用於排序矩形
bool compareRect(const Rect &a, const Rect &b) {
    if (a.yStart != b.yStart)
        return a.yStart < b.yStart; // 先按照 y 座標排序（從下到上）
    else
        return a.xStart < b.xStart; // 若 y 相同，則按照 x 座標排序（從左到右）
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
    int numnodes = movableNodes + numTerminals;

    readNetsFile(benchmarkPath, netsFile);
    readPlFile(benchmarkPath, plFile);
    readSclFile(benchmarkPath, sclFile);

    /* -----------------------------------------------------------------------------
        Modified Tetris Legalization Algorithm
    -------------------------------------------------------------------------------- */

    // 定義新的座標向量
    std::vector<float> xCellCoord_p(numnodes + 1);
    std::vector<float> yCellCoord_p(numnodes + 1);

    // 構建矩形列表
    std::vector<Rect> rects;
    for (int i = 1; i <= numnodes; ++i) {
        Rect rect;
        rect.xStart = xCellCoord[i];
        rect.yStart = yCellCoord[i];
        rect.xEnd = rect.xStart + cellWidth[i];
        rect.yEnd = rect.yStart + cellHeight[i];
        rect.index = i;
        rects.push_back(rect);
    }

    // 按照 y 座標從下到上，x 座標從左到右排序
    std::sort(rects.begin(), rects.end(), compareRect);

    // 建立每個 row 的可用區間列表
    std::vector<std::vector<std::pair<float, float>>> availableIntervals(numRows);
    for (int row = 0; row < numRows; ++row) {
        availableIntervals[row].push_back({rowOriginX[row], rowEndX[row]});
    }

    // 開始 Modified Tetris Legalization
    for (auto &rect : rects) {
        float width = cellWidth[rect.index];
        float height = cellHeight[rect.index];

        // 儲存最佳位置資訊
        bool placed = false;
        float minDistance = std::numeric_limits<float>::max();
        int bestRow = -1;
        float bestXStart = rect.xStart;
        float bestYStart = rect.yStart;

        // 遍歷所有 row，尋找最適合的位置
        for (int row = 0; row < numRows; ++row) {
            float yStart = siteOriginY + row * coreRowHeight;

            // 確保矩形高度不超過 row 的高度
            if (height > coreRowHeight) {
                continue; // 無法放置在該 row
            }

            // 在該 row 中遍歷可用區間
            auto &intervals = availableIntervals[row];
            for (auto it = intervals.begin(); it != intervals.end(); ++it) {
                float intervalStart = it->first;
                float intervalEnd = it->second;

                // 如果區間足夠放下矩形
                if (intervalEnd - intervalStart >= width) {
                    // 計算放置位置，盡量接近原始 x 和 y 座標
                    float xStart = clamp(rect.xStart, intervalStart, intervalEnd - width);
                    float distance = std::abs(xStart - rect.xStart) + std::abs(yStart - rect.yStart);

                    if (distance < minDistance) {
                        minDistance = distance;
                        bestXStart = xStart;
                        bestYStart = yStart;
                        bestRow = row;
                        placed = true;
                    }
                }
            }
        }

        if (placed) {
            // 更新最佳 row 的可用區間
            float xPlacedStart = bestXStart;
            float xPlacedEnd = xPlacedStart + width;

            auto &intervals = availableIntervals[bestRow];
            std::vector<std::pair<float, float>> newIntervals;

            for (auto &interval : intervals) {
                if (interval.second <= xPlacedStart || interval.first >= xPlacedEnd) {
                    // 不重疊，保留
                    newIntervals.push_back(interval);
                } else {
                    // 可能需要拆分區間
                    if (interval.first < xPlacedStart) {
                        newIntervals.push_back({interval.first, xPlacedStart});
                    }
                    if (interval.second > xPlacedEnd) {
                        newIntervals.push_back({xPlacedEnd, interval.second});
                    }
                }
            }
            intervals = newIntervals;

            // 記錄新位置
            xCellCoord_p[rect.index] = bestXStart;
            yCellCoord_p[rect.index] = bestYStart;
        } else {
            // 無法放置，保持原始位置並輸出警告
            xCellCoord_p[rect.index] = rect.xStart;
            yCellCoord_p[rect.index] = rect.yStart;
            std::cout << "Warning: Rect " << rect.index << " could not be placed." << std::endl;
        }
    }

    // 輸出結果
    // for (int i = 1; i <= numnodes; ++i) {
    //     std::cout << "Rect " << i << " new position: (" << xCellCoord_p[i] << ", " << yCellCoord_p[i] << ")" << std::endl;
    // }

    // 如果需要，將結果寫入檔案
    write_python_File( xCellCoord_p.data(), yCellCoord_p.data());


    //計算總共cost以及最大的cost
    float total_cost=0.0;
    float max_cost=0.0;
    for(int i=1;i<xCellCoord_p.size();i++){
        float ab_x=xCellCoord_p[i]-xCellCoord[i];
        float ab_y=yCellCoord_p[i]-yCellCoord[i];
        // printf("old xCellCoord[i]: %.2f yCellCoord[i]: %.2f\n",xCellCoord[i],yCellCoord[i]);
        // printf("new xCellCoord_p[i]: %.2f yCellCoord_p[i]: %.2f\n",xCellCoord_p[i],yCellCoord_p[i]);
        float cost= std::fabs(xCellCoord_p[i]-xCellCoord[i])+std::fabs(yCellCoord_p[i]-yCellCoord[i]);
        total_cost+=cost;
        if(max_cost<cost){
            max_cost=cost;
        }
    }
    printf("total_cost: %.2f\n",total_cost);
    printf("max_cost  : %.2f\n",max_cost);

    // 檢查是否有重疊的矩形
    bool hasOverlap = false;
    for (int i = 1; i <= numnodes; ++i) {
        for (int j = i + 1; j <= numnodes; ++j) {
            // 獲取矩形 i 的範圍
            float xStart_i = xCellCoord_p[i];
            float xEnd_i = xStart_i + cellWidth[i];
            float yStart_i = yCellCoord_p[i];
            float yEnd_i = yStart_i + cellHeight[i];

            // 獲取矩形 j 的範圍
            float xStart_j = xCellCoord_p[j];
            float xEnd_j = xStart_j + cellWidth[j];
            float yStart_j = yCellCoord_p[j];
            float yEnd_j = yStart_j + cellHeight[j];

            // 檢查 X 軸上的重疊，邊界相等不算重疊
            bool xOverlap = (xStart_i < xEnd_j && xEnd_i > xStart_j)||
                            (xStart_i > xEnd_j && xEnd_i < xStart_j);
            // 如果邊界相等，不算重疊
            if (xEnd_i == xStart_j || xEnd_j == xStart_i) {
                xOverlap = false;
            }

            // 檢查 Y 軸上的重疊，邊界相等不算重疊
            bool yOverlap = yStart_i < yEnd_j && yEnd_i > yStart_j;
            // if (yEnd_i == yStart_j || yEnd_j == yStart_i) {
            //     yOverlap = false;
            // }

            if (xOverlap && yOverlap) {
                hasOverlap = true;
                std::cout << "Overlap detected between Rect " << i <<" ("<<xStart_i<<" , "<<yStart_i<< " ) and Rect " 
                          << j <<" ("<<xStart_j<<" , "<<yStart_j<< " )" << std::endl;
            }
        }
    }

    if (!hasOverlap) {
        std::cout << "No overlaps detected." << std::endl;
    }


    freeHash();
    return 0;
}
