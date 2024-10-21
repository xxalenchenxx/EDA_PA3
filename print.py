import matplotlib.pyplot as plt

# 讀取資料並將矩形繪製出來
def plot_rectangles(file_path):
    # 打開檔案並讀取資料
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # 讀取範圍
    siteOriginX, siteOriginY, siteEndX, siteEndY = map(float, lines[0].split())

    # 創建圖表
    fig, ax = plt.subplots()

    # 迭代每一行資料並繪製矩形外框 (從第二行開始讀取)
    for line in lines[1:]:
        # 解析每行資料
        x, y, width, height = map(float, line.split())

        # 繪製矩形 (只畫外框)
        rect = plt.Rectangle((x, y), width, height, fill=False, edgecolor='blue', linewidth=2)
        ax.add_patch(rect)

    # 設定x, y軸的標籤和範圍
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_aspect('equal', adjustable='box')

    # 設定顯示範圍為整張圖的範圍
    ax.set_xlim([siteOriginX, siteEndX])
    ax.set_ylim([siteOriginY, siteEndY])

    # 保存圖片為 output.jpg
    plt.savefig('output.jpg', format='jpg')

    # 如果你還是想要顯示圖片，可以選擇性地保留這一行
    # plt.show()

# 設定檔案路徑
file_path = 'plt.txt'

# 呼叫函數來繪製並保存矩形外框
plot_rectangles(file_path)
