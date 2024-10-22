# 編譯器與優化選項
CC = g++ -O

# 可執行檔案名稱
TARGET = PA3

# 需要編譯的檔案
SRCS = PA3.cpp memAlloc.c bookshelf_IO.c

# 物件檔案
OBJS = PA3.o memAlloc.o bookshelf_IO.o

# 編譯可執行檔
all: $(TARGET)

# 主程式及依賴檔案的編譯和鏈接
$(TARGET): $(OBJS)
	$(CC) -o $(TARGET) $(OBJS) -lm

# 各個檔案的編譯規則
PA3.o: PA3.cpp
	$(CC) -c PA3.cpp -o PA3.o

memAlloc.o: memAlloc.c memAlloc.h
	$(CC) -c memAlloc.c -o memAlloc.o

bookshelf_IO.o: bookshelf_IO.c bookshelf_IO.h
	$(CC) -c bookshelf_IO.c -o bookshelf_IO.o

# 清理物件檔和可執行檔
clean:
	rm -f *.o $(TARGET) output.jpg
