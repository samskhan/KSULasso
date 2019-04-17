from RandomLasso import RandomLasso
from sklearn import linear_model



def main():
    rl = RandomLasso(linear_model.Lasso(alpha=0),"/home/delilah2/Documents/KSULasso/KSULasso_Python/sample/test2_x.csv","/home/delilah2/Documents/KSULasso/KSULasso_Python/sample/test2_y.csv")
    #rl.readData("/home/delilah2/Documents/KSULasso/KSULasso_Python/sample/test2_x.csv","/home/delilah2/Documents/KSULasso/KSULasso_Python/sample/test2_y.csv")
    rl.printCoeff()
    
if __name__== "__main__":
    main()