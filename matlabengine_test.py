import matlab.engine
eng = matlab.engine.start_matlab()
s = eng.genpath('C:/downloads_/lab/lab_taulu')
eng.addpath(s, nargout=0)

# attempt to load in a matrix from matlab
test_matrix = eng.test_matlab_function_foster()
print(test_matrix)
