import cv2
import time
import copy

for i in range(1, 4):
    print('loading upper jaw number %d' % i)
    img_address = './upper_jaws/%d_upper_clahe_sauvola.bmp' % i
    img = cv2.imread(img_address, 0)
    img_copy = copy.deepcopy(x=img)
    
    print('applying genetic algorithm inorder to find best tooth separator lines')
    t0 = time.time()
    # draw tooth separator lines 
    t1 = time.time()
    print('elapsed time for GA: %.2f secs' % (t1 - t0))

    print('saving result image')
    # cv2.imwrite('./upper_jaws/%d_upper_clahe_sauvola.bmp' % i, final_image)
    print('------------------------ result saved ------------------------')
    

for i in range(1, 4):
    print('loading lower jaw number %d' % i)
    img_address = './lower_jaws/%d_lower_clahe_sauvola.bmp' % i
    img = cv2.imread(img_address, 0)
    img_copy = copy.deepcopy(x=img)
    
    print('applying genetic algorithm inorder to find best tooth separator lines')
    t0 = time.time()
    # draw tooth separator lines 
    t1 = time.time()
    print('elapsed time for GA: %.2f secs' % (t1 - t0))

    print('saving result image')
    # cv2.imwrite('./lower_jaws/%d_lower_clahe_sauvola.bmp' % i, final_image)
    print('------------------------ result saved ------------------------')


print('process finished')