from sklearn import svm

def fit_and_score(split_exs): #[ex_train, ex_test]
    exs_train,exs_test = split_exs
    classer = fit_svm(exs_train)
    tested = score_svm(classer, exs_test)
    return tested

def fit_svm(exs_train):
    classer = svm.SVC(kernel='linear', probability=True) 
    _,_,X,y = extract_data_labels(exs_train.examples)
    classer.fit(X,y)
    return classer

def score_svm(classer, exs_test):
    p1s,p2s,X0,y0 = extract_data_labels(exs_test.examples)
    y0p = [x[1] for x in classer.predict_proba(X0)]
    tested = [(p1,p2,score,label) for p1,p2,score,label in zip(p1s,p2s,y0p,y0)]
    tested.sort(key=lambda x:x[2],reverse=True)
    return tested

def extract_data_labels(ex_list):
    p1 = [e[0] for e in ex_list]
    p2 = [e[1] for e in ex_list]
    X = [[0 if x == '?' else x for x in e[3:]] for e in ex_list]
    y = [int(e[2]=='true') for e in ex_list]
    return p1,p2,X,y

def rescale(p,n):
    """
    Rescale posterior probability p according to multiple of negatives n.
    """
    return p/(1+(1-p)*n)
