import numpy as np

def crude_match(array1,array2,criteria,verbose=False):
    '''
    Conduct crude match between two given list of positions of stars
    
    INPUTS:
    array1: the first star array with first and second column array1[:,0:2] are star position(x,y) and array1.shape = (N1,m)
    array2: the second star array with first and second column array1[:,0:2] are star position(x,y)and array1.shape = (N2,m)
            array1 and array have the same width(columns) and ONE of N1 and N2 can equal 1
    criteria: = c(PIXES) the criteria for match is that one star is within c pixes from the corrresponding star's central position
    '''
    
    if len(array1) < len(array2):
        refarray = array2
        matarray = array1
        alter = 1
    else:
        refarray = array1
        matarray = array2
        alter = 0
    
    N1 = len(matarray)
    N2 = len(refarray)
        
    matmask=[]
    refmask=[]
    
    if matarray.shape == (2,) or matarray.shape == (1,2):
	matarray = matarray.copy().ravel()
	#print matarray.shape
	#print type(refarray)
	#print refarray.shape
        diffarray = refarray - matarray        
        temp = []
        for j,diff in enumerate(diffarray):
            #print diff
            #print diff[0]		
	    #print diff
            #print diff.shape
    	    #print diff[0]
	    #print diff[1]
            if np.abs(diff[0])<criteria and np.abs(diff[1])< criteria:                
                temp.append(j)
                #print diff            
        
        if len(temp)==1:
            matched_ref = refarray[temp[0]]            
            matched_mat = matarray
            matmask.append(1)
            refmask.append(temp[0])
        else:
            if len(temp)>1:
                print 'more than one objects fall into the criteria region of reference star. This object is droped!'
            stars1 = []
            stars2 = []
            mask1 =[]
            mask2 =[]
            return stars1,mask1,stars2,mask2            
    else:        
        for i,current in enumerate(matarray):                    
            diffarray = refarray - current
            if verbose:        
                print diffarray
            temp = []
            for j,diff in enumerate(diffarray):
                try:
                    diff1 = diff[0]
                    diff2 = diff[1]
                except:
                    diff1 = diff[0,0]
                    diff2 = diff[0,1]
                    
                if np.abs(diff1)<criteria and np.abs(diff2)< criteria:
                    temp.append(j)
                    #print diff
                
            if len(temp)==1:
                matmask.append(i)
                refmask.append(temp[0])
            elif len(temp)>1:
                print 'more than one objects fall into the criteria region of reference star. This object is droped!'
    
        if verbose:        
            print refmask
            print matmask
        
        matched_ref = refarray[refmask]
        matched_mat = matarray[matmask]
        
    #make sure the output stars1,stars2 correspond to array1,array2 respectively
    
    if alter:
        stars1 = matched_mat
        stars2 = matched_ref
        mask1 = matmask
        mask2 = refmask
    else:
        stars1 = matched_ref
        stars2 = matched_mat
        mask1 = refmask
        mask2 = matmask    
    
    if verbose:
        print "From array1:\n"
        print stars1
        print "From array2:\n"
        print stars2
    
    return stars1,mask1,stars2,mask2


def find_corresponding(basearray,target,criteria,verbose=False):
    '''
	find the location of 'target' in 'array' when it meet the 'criteria'
    '''
    target = target.reshape((1,2))
    diffs = basearray-target
    ret_mask= []
    for i,diff in enumerate(diffs):
        if np.abs(diff[0])<criteria and np.abs(diff[1])<criteria:
	    ret_mask.append(i)
    if ret_mask == []:
	exist = False
	return exist,0,0
    else:
	exist = True 
        return exist,basearray[ret_mask],ret_mask
    
