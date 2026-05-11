#ifndef REMAPKNN_H
#define REMAPKNN_H

#include <string>
#include "knndata.h"

struct RemapknnParams
{
  std::string gridString;
  KnnParams knnParams;
};

RemapknnParams remapknn_get_parameter();
void remapknn_verify_parameter(KnnParams const &knnParams);
void print_knn_parameter(KnnParams const &knnParams, std::string const &prefix);
void remapknn_print_parameter(RemapknnParams const &params);

#endif
