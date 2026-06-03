import numpy as np

def forward_difference(data: np.ndarray, i: int):
    return data[i+1] - data[i]

def symmetric_difference(data: np.ndarray, i: int):
    return (data[i+1] - data[i-1]) / 2

def derivative(data: np.ndarray):
    derivative = np.zeros_like(data)
    for i in range(1, len(data) - 1):
        derivative[i] = symmetric_difference(data, i)

    derivative[len(data) - 1] = forward_difference(data, len(data) - 2)

    return derivative