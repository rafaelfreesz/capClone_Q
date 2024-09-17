//
// Created by rafael on 22/07/22.
//

#ifndef CAPCLONE_UTILS_H
#define CAPCLONE_UTILS_H


#include "Instance.h"

class Utils {
public:
    static void tokenize(string str, vector<string> &token_v, string DELIMITER);
    static vector<Instance *> loadInstances(string strFileName);
    static void testInstance(vector<Instance*> instances);
    static void swapInt(int i, int j, int* array);
};


#endif //CAPCLONE_UTILS_H
