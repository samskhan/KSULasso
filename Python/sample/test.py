from RandomLasso import RandomLasso
from sklearn import linear_model



def main():

    #rl = RandomLasso(linear_model.Lasso(alpha=0),"/home/delilah2/Documents/KSULasso/KSULasso_Python/sample/test2_x.csv","/home/delilah2/Documents/KSULasso/KSULasso_Python/sample/test2_y.csv")
    #rl.readData("/home/delilah2/Documents/KSULasso/KSULasso_Python/sample/test2_x.csv","/home/delilah2/Documents/KSULasso/KSULasso_Python/sample/test2_y.csv")
    rl = RandomLasso()
    x,y=rl.readData("/home/delilah2/Documents/KSULasso/KSULasso_Python/sample/test2_x.csv","/home/delilah2/Documents/KSULasso/KSULasso_Python/sample/test2_y.csv")
    # rl.printData()
    print(rl.bootstraps)
    print(rl.sampleCount)
    print(rl.featureCount)

    coeff = rl.fit(linear_model.Lasso(alpha=1.0),linear_model.Lasso(alpha=1.0),x=x,y=y)
    
    print(coeff)
if __name__== "__main__":
    main()