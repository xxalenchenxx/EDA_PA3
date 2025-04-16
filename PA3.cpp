#include <string.h>
// #include <stdio.h>
// #include <stdlib.h>
#include <map>
#include <cmath>
#include <ctime>
#include <vector>
#include <chrono> 
#include <limits>
#include <cstdlib>
#include <iostream>
#include <algorithm>        
#include "bookshelf_IO.h"

double time_limit=120.0; //set limit time of SA tiime

struct Rect {
    float xStart, xEnd;
    float yStart, yEnd;
    int index; 
};



float clamp(float val, float minVal, float maxVal) {
    return std::max(minVal, std::min(maxVal, val));
}


bool compareRect(const Rect &a, const Rect &b) {
    float sumA = (a.xEnd - a.xStart) + (a.yEnd - a.yStart);
    float sumB = (b.xEnd - b.xStart) + (b.yEnd - b.yStart);
    return sumA > sumB;
}

void check_overlap(std::vector<float> xCellCoord_p,std::vector<float> yCellCoord_p){
    
    bool hasOverlap = false;
    for (int i = 1; i <= numNodes; ++i) {
        for (int j = i + 1; j <= numNodes; ++j) {
            
            float xStart_i = xCellCoord_p[i];
            float xEnd_i = xStart_i + cellWidth[i];
            float yStart_i = yCellCoord_p[i];
            float yEnd_i = yStart_i + cellHeight[i];

            
            float xStart_j = xCellCoord_p[j];
            float xEnd_j = xStart_j + cellWidth[j];
            float yStart_j = yCellCoord_p[j];
            float yEnd_j = yStart_j + cellHeight[j];

            
            bool xOverlap = xStart_i < xEnd_j && xEnd_i > xStart_j;
            if (xEnd_i == xStart_j || xEnd_j == xStart_i) {
                xOverlap = false;
            }

           
            bool yOverlap = yStart_i < yEnd_j && yEnd_i > yStart_j;
            if (yEnd_i == yStart_j || yEnd_j == yStart_i) {
                yOverlap = false;
            }

            if (xOverlap && yOverlap) {
                hasOverlap = true;
                std::cout << "Overlap detected between Rect " << cellName[i] <<" ("<<xStart_i<<" , "<<yStart_i<< " ) and Rect " 
                          << cellName[j] <<" ("<<xStart_j<<" , "<<yStart_j<< " )" << std::endl;
            }
        }
    }

    if (hasOverlap) {
        std::cout << "[ERROR]have overlaps" << std::endl;
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
        // 檢查 row 是否合法
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
    if (has_non_site) {
        std::cout << "[ERROR]positions are not on site grid" << std::endl;
    }

}


bool checkOverlap(const std::vector<float> &xCellCoord_p, const std::vector<float> &yCellCoord_p, int numNodes) {
    for (int i = 1; i <= numNodes; ++i) {
        for (int j = i + 1; j <= numNodes; ++j) {
            
            float xStart_i = xCellCoord_p[i];
            float xEnd_i = xStart_i + cellWidth[i];
            float yStart_i = yCellCoord_p[i];
            float yEnd_i = yStart_i + cellHeight[i];

            
            float xStart_j = xCellCoord_p[j];
            float xEnd_j = xStart_j + cellWidth[j];
            float yStart_j = yCellCoord_p[j];
            float yEnd_j = yStart_j + cellHeight[j];

            
            bool xOverlap = xStart_i < xEnd_j && xEnd_i > xStart_j;
            if (xEnd_i == xStart_j || xEnd_j == xStart_i) {
                xOverlap = false;
            }

            
            bool yOverlap = yStart_i < yEnd_j && yEnd_i > yStart_j;
            if (yEnd_i == yStart_j || yEnd_j == yStart_i) {
                yOverlap = false;
            }

            if (xOverlap && yOverlap) {
                return true; 
            }
        }
    }
    return false; 
}


void calculateDisplacement(const std::vector<float> &xCellCoord_p, const std::vector<float> &yCellCoord_p, int numNodes, float &total_cost) {
    total_cost = 0.0f;
    float max_cost = 0.0f;
    for (int i = 1; i <= numNodes; ++i) {
        float ab_x = xCellCoord_p[i] - xCellCoord[i];
        float ab_y = yCellCoord_p[i] - yCellCoord[i];
        float cost = std::fabs(ab_x) + std::fabs(ab_y);
        total_cost += cost;
        if(cost>max_cost){
            max_cost=cost;
        }
    }
    total_cost+=max_cost;
}


void simulatedAnnealing(std::vector<float> &xCellCoord_p, std::vector<float> &yCellCoord_p,
    float *cellWidth, float *cellHeight,
    const std::vector<float> &siteSpacing, const std::vector<float> &subrowOrigin,
    int numNodes, int numRows, float siteOriginY, float coreRowHeight,auto start) {

    float T = 1024.0f; // intial T
    float T_min = 1e-1f; // minimal T
    float alpha = 0.995f; // coe
    int maxIterations = 2048; // max iteration


    std::map<float, std::vector<int>> widthToCells;
    for (int i = 1; i <= numNodes; ++i) {
        widthToCells[cellWidth[i]].push_back(i);
    }

    std::vector<float> currentX = xCellCoord_p;
    std::vector<float> currentY = yCellCoord_p;


    float currentTotalCost, currentMaxCost;
    calculateDisplacement(currentX, currentY, numNodes, currentTotalCost);
    float historyTotalCost=currentTotalCost;
    std::srand(static_cast<unsigned int>(std::time(nullptr))); // random
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    //start SA
    while (T > T_min && elapsed.count() < time_limit) {
    // std::cout<<"T:"<<T<<std::endl;
    // std::cout<<"end: "<<elapsed.count()<<std::endl;
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
        end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;
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
    auto start = std::chrono::high_resolution_clock::now();
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
        // if(i<10)
        //     printf("%s\t%f\t%f\n",cellName[i],xCellCoord[i],yCellCoord[i]);
    }
    // write_python_File( xCellCoord, yCellCoord);

    std::sort(rects.begin(), rects.end(), compareRect);

    // for(auto &a:rects){
    //     printf("a.x: %f a.y: %f\n",a.xStart,a.yStart);
    // }

    std::vector<std::vector<std::pair<float, float>>> availableIntervals(numRows);
    std::vector<std::vector<std::pair<float, float>>> availableIntervals_base(numRows);
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
    availableIntervals_base=availableIntervals;
    bool base_flag=false;
    // Modified Tetris Legalization
    for (auto &rect : rects) {
        // 假設已在迴圈中取得目前 cell 的寬度與高度：
        float width = cellWidth[rect.index];
        float height = cellHeight[rect.index];

        bool placed = false;
        float minDistance = std::numeric_limits<float>::max();
        int bestRow = -1;
        float bestXStart = rect.xStart;
        float bestYStart = rect.yStart;

        for (int row = 0; row < numRows; ++row) {
            float yStart = siteOriginY + row * coreRowHeight;

            // 若 cell 太高無法放入此 row 則跳過
            if (height > coreRowHeight) {
                continue;
            }

            float spacing = siteSpacing[row];
            float origin = subrowOrigin[row]; // 該 row 的起始 x 座標
            int sites = numSites[row];         // 該 row 可用的站點數量

            auto &intervals = availableIntervals[row];

            // 若該 row 只有15個可用區段：使用原本逐站點掃描搜尋最佳候選位置
            if (intervals.size() < 15) {
                for (auto it = intervals.begin(); it != intervals.end(); ++it) {
                    float intervalStart = it->first;
                    float intervalEnd = it->second;

                    int startSiteIndex = std::max(0, (int)std::ceil((intervalStart - origin) / spacing));
                    int endSiteIndex = std::min(sites, (int)((intervalEnd - origin) / spacing));

                    for (int siteIndex = startSiteIndex; siteIndex < endSiteIndex; ++siteIndex) {
                        float xSite = origin + siteIndex * spacing;
                        // 若該候選位置放置 cell 超出該間隔則中斷
                        if (xSite + width >= intervalEnd) {
                            break;
                        }
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
            // 若該 row 有多個區段：每個區段只評估最左側與最右側候選位置
            else {
                for (auto it = intervals.begin(); it != intervals.end(); ++it) {
                    float intervalStart = it->first;
                    float intervalEnd = it->second;

                    int startSiteIndex = std::max(0, (int)std::ceil((intervalStart - origin) / spacing));
                    int endSiteIndex = std::min(sites, (int)((intervalEnd - origin) / spacing));

                    // 評估左邊界候選位置
                    float candidateLeft = origin + startSiteIndex * spacing;
                    if (candidateLeft + width < intervalEnd) {
                        float distance = std::fabs(candidateLeft - rect.xStart) + std::fabs(yStart - rect.yStart);
                        if (distance < minDistance) {
                            minDistance = distance;
                            bestXStart = candidateLeft;
                            bestYStart = yStart;
                            bestRow = row;
                            placed = true;
                        }
                    }

                    // 評估右邊界候選位置 (注意：endSiteIndex 為區段中不包含的上界，故使用 endSiteIndex - 1)
                    if (endSiteIndex > startSiteIndex) {
                        float candidateRight = origin + (endSiteIndex - 1) * spacing;
                        if (candidateRight + width < intervalEnd) {
                            float distance = std::fabs(candidateRight - rect.xStart) + std::fabs(yStart - rect.yStart);
                            if (distance < minDistance) {
                                minDistance = distance;
                                bestXStart = candidateRight;
                                bestYStart = yStart;
                                bestRow = row;
                                placed = true;
                            }
                        }
                    }
                }
            }
        }

        if (placed) {
            // 已找到合適位置，更新該 row 的可用區間
            float xPlacedStart = bestXStart;
            float xPlacedEnd = xPlacedStart + width;

            auto &intervals = availableIntervals[bestRow];
            std::vector<std::pair<float, float>> newIntervals;
            for (auto &interval : intervals) {
                if (interval.second <= xPlacedStart || interval.first >= xPlacedEnd) {
                    newIntervals.push_back(interval);
                } else {
                    if (interval.first < xPlacedStart) {
                        newIntervals.push_back({interval.first, xPlacedStart});
                    }
                    if (interval.second > xPlacedEnd) {
                        newIntervals.push_back({xPlacedEnd, interval.second});
                    }
                }
            }
            intervals = newIntervals;

            xCellCoord_p[rect.index] = bestXStart;
            yCellCoord_p[rect.index] = bestYStart;
        } else {
            // 若無法找到合法位置，則保持原位並發出警告訊息
            xCellCoord_p[rect.index] = rect.xStart;
            yCellCoord_p[rect.index] = rect.yStart;
            std::cout << "Warning: Rect " << cellName[rect.index] 
                      << " could not be placed. (" 
                      << xCellCoord_p[rect.index] << " , " 
                      << yCellCoord_p[rect.index] << " )" 
                      << std::endl;
            base_flag = true;
            break;
        }

    }


    float total_cost=0.0;
    float max_cost=0.0;
    // for(int i=1;i<xCellCoord_p.size();i++){
    //     float ab_x=xCellCoord_p[i]-xCellCoord[i];
    //     float ab_y=yCellCoord_p[i]-yCellCoord[i];
    //     float cost= std::fabs(ab_x)+std::fabs(ab_y);
    //     total_cost+=cost;
    //     if(max_cost<cost){
    //         max_cost=cost;
    //     }
    // }
    // printf("[before]Total displacement: %.1f\n",total_cost);
    // printf("[before]Maximum displacement: %.1f\n",max_cost);

    /* -----------------------------------------------------------------------------
        simulated Annealing Algorithm 
    -------------------------------------------------------------------------------- */
    if(base_flag){
        // 再次嘗試放置未放置的矩形，逐行進行緊密排列
       // 對於每一個矩形，根據 compareRect 的順序重新從下到上、從左到右進行排放
       printf("do base case\n");
        for (auto &rect : rects) {
            for (int row = 0; row < numRows; ++row) {
                float yStart = siteOriginY + row * coreRowHeight;

                // 如果矩形的高度超過當前 row 的高度，則跳過
                // if ((rect.yEnd - rect.yStart) > coreRowHeight) {
                //     continue;
                // }

                // 計算該行的間隔參數
                float spacing = siteSpacing[row];
                float origin = subrowOrigin[row];
                int sites = numSites[row];
                auto &intervals = availableIntervals_base[row];

                // 嘗試在當前行的可用間隔中找到合適的區域來放置矩形
                for (auto it = intervals.begin(); it != intervals.end(); ++it) {
                    float intervalStart = it->first;
                    float intervalEnd = it->second;

                    // 計算可以開始放置的位置的索引
                    int startSiteIndex = std::max(0, (int)std::ceil((intervalStart - origin) / spacing));
                    int endSiteIndex = std::min(sites, (int)((intervalEnd - origin) / spacing));

                    // 在當前可用的 site 中進行檢查
                    for (int siteIndex = startSiteIndex; siteIndex < endSiteIndex; ++siteIndex) {
                        float xSite = origin + siteIndex * spacing;

                        // 確保矩形可以完整地放置在這個區間內
                        if (xSite + (rect.xEnd - rect.xStart) > intervalEnd) {
                            break; // 無法在該區間放置該矩形
                        }

                        // 記錄最佳放置位置
                        float xPlacedStart = xSite;
                        float xPlacedEnd = xPlacedStart + (rect.xEnd - rect.xStart);

                        // 更新 availableIntervals，移除已經放置的位置
                        std::vector<std::pair<float, float>> newIntervals;
                        for (auto &interval : intervals) {
                            if (interval.second <= xPlacedStart || interval.first >= xPlacedEnd) {
                                // 如果區間在放置位置的左邊或右邊，則保留
                                newIntervals.push_back(interval);
                            } else {
                                // 如果區間被部分覆蓋，則進行劃分
                                if (interval.first < xPlacedStart) {
                                    newIntervals.push_back({interval.first, xPlacedStart});
                                }
                                if (interval.second > xPlacedEnd) {
                                    newIntervals.push_back({xPlacedEnd, interval.second});
                                }
                            }
                        }
                        intervals = newIntervals;

                        // 更新矩形的放置位置
                        xCellCoord_p[rect.index] = xPlacedStart;
                        yCellCoord_p[rect.index] = yStart;

                        // 矩形已經成功放置，跳出當前矩形的處理
                        goto next_rect;
                    }
                }
            }

            // 如果矩形未能被放置，給出警告訊息
            std::cout << "Warning: Rect " << cellName[rect.index] << " could not be placed after reordering."
                      << " (" << rect.xStart << ", " << rect.yStart << ")" << std::endl;

        next_rect:;
        }
    
    }else{
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        while(elapsed.count()<time_limit){
            // printf("first: %d\n",first);
            simulatedAnnealing(xCellCoord_p, yCellCoord_p, cellWidth, cellHeight, siteSpacing, subrowOrigin,
                               numNodes, numRows, siteOriginY, coreRowHeight,start);

            end = std::chrono::high_resolution_clock::now();
            elapsed = end - start;
        }
    }
    




    // visualize result if needed
    // write_python_File( xCellCoord_p.data(), yCellCoord_p.data());


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
    printf("\n-------------result-------------\n");
    printf("Total displacement: %.1f\n",total_cost);
    printf("Maximum displacement: %.1f\n",max_cost);

    // check_overlap
    check_overlap(xCellCoord_p,yCellCoord_p);
    //檢查是否符合 site的座標
    check_on_site(xCellCoord_p,yCellCoord_p,siteSpacing,subrowOrigin);

    //output file
    clear_output_dir(outputPath);
    write_aux_File(outputPath,benchmarkName);
    writePlFile(outputPath,benchmarkName,xCellCoord_p.data(),yCellCoord_p.data());
    write_other_File(benchmarkPath,outputPath,benchmarkName);

    freeHash();
    return 0;
}
