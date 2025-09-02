# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 14:34:30 2023

@author: Cheng-Hao Yu
"""

#%% 第一題：計算兩個三維向量的和
# 給定兩個向量，請寫一個函數 add_vectors(vector1, vector2) 來計算它們的和。

def add_vectors(vector1, vector2):
    result = [vector1[i] + vector2[i] for i in range(len(vector1))]
    return result

# Example Usage
vector1 = [2, 2, 1]
vector2 = [3, 3, 2]
result1 = add_vectors(vector1, vector2)


#%% 第二題：計算兩個三維向量的外積
# 請撰寫一個名為 vector_outer_product(vector1, vector2) 的函式，接受兩個三維向量作為輸入，並返回它們的外積。

def vector_outer_product(vector1, vector2):
    result = [0, 0, 0]
    result[0] = vector1[1]*vector2[2] - vector1[2]*vector2[1]
    result[1] = vector1[2]*vector2[0] - vector1[0]*vector2[2]
    result[2] = vector1[0]*vector2[1] - vector1[1]*vector2[0]
    return result

# Example Usage
vector1 = [2, 2, 1]
vector2 = [3, 3, 2]
result2 = vector_outer_product(vector1, vector2)


#%% 第三題：計算轉置矩陣
# 給定一個矩陣 matrix，請寫一個函式 transpose(matrix) 將其進行轉置。轉置後的矩陣應該是將原矩陣的行變成列，列變成行。

def transpose(matrix):
    rows = len(matrix)
    cols = len(matrix[0])
    
    transposed = []
    for i in range(cols):
        transposed.append([0] * rows)

    for j in range(rows):
        for k in range(cols):
            transposed[k][j] = matrix[j][k]
    
    return transposed

# Example Usage
matrix1 = [[1, 2, 3],[4, 5, 6],[7, 8, 9]]
result3 = transpose(matrix1)


#%% 第四題：計算兩個向量的內積
# 請寫一個 Python 函數 vector_dot_product(vector1, vector2)，接受兩個等長的向量，計算它們的內積。若兩個列表的長度不相等，請返回 None。

def vector_dot_product(vector1, vector2):
    if len(vector1) != len(vector2):
        return None

    result = 0

    for i in range(len(vector1)):
        result += vector1[i] * vector2[i]

    return result

# Example Usage
vector1 = [2, 2, 1]
vector2 = [3, 3, 2]
result4 = vector_dot_product(vector1, vector2)


#%% 第五題：計算單位向量
# 請寫一個 Python 函式 normalize_vector(vector)，接受一個包含整數的向量 vector，並回傳該向量的單位向量。向量的長度為 0 時無法進行單位向量的計算，此時請回傳 None。

def normalize_vector(vector):
    length_squared = sum(x ** 2 for x in vector)
    
    if length_squared == 0:
        return None
    
    length = length_squared**0.5
    unit_vector = [x / length for x in vector]
    return unit_vector

# Example Usage
vector1 = [2, 2, 1]
result5 = normalize_vector(vector1)