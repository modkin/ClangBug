#include <iostream>
#include <cmath>

inline size_t num_microvertices_per_edge(size_t level)
{
  return (size_t) std::pow(2, level) + 1;
}


inline size_t num_microvertices_per_face(size_t level)
{
  return (size_t) ((std::pow(2, level)+1) * (std::pow(2, level-1) + 1));
}

inline void prolongate(size_t level)
{
  size_t rowsize_coarse = num_microvertices_per_edge(level);
  size_t rowsize_fine = num_microvertices_per_edge(level+1);

  size_t totalPoints =  num_microvertices_per_face(level);


  double* face_data_f = (double *) malloc(num_microvertices_per_face(level-2+1)*sizeof(double));
  for(int i = 0; i < num_microvertices_per_face(level-2+1)*sizeof(double); ++i){
    face_data_f[i] = 0.1;
  }
  double* face_data_c = (double *) malloc(num_microvertices_per_face(level-2)*sizeof(double));
  for(int i = 0; i < num_microvertices_per_face(level-2)*sizeof(double); ++i){
    face_data_c[i] = 5;
  }


  size_t mr_c = 1;
  size_t mr_f = rowsize_fine + 2;

  size_t i_rowsize_coarse = rowsize_coarse;

  for (size_t i_coarse = 0; i_coarse < rowsize_coarse-2; ++i_coarse)
  {
    for (size_t j_coarse = 0; j_coarse < i_rowsize_coarse-3; ++j_coarse)
    {
      face_data_f[mr_f] = 0.5 * (face_data_c[mr_c] + face_data_c[mr_c + i_rowsize_coarse]);
      face_data_f[mr_f-1] = 0.5 * (face_data_c[mr_c] + face_data_c[mr_c + i_rowsize_coarse - 1]);
      //size_t OMG = mr_f + rowsize_fine - 1 - 1;
      face_data_f[mr_f + rowsize_fine - 1 - 1] = 0.5 * (face_data_c[mr_c + i_rowsize_coarse] + face_data_c[mr_c + i_rowsize_coarse - 1]);

      face_data_f[mr_f + rowsize_fine - 1] = face_data_c[mr_c + i_rowsize_coarse];

      mr_c += 1;
      mr_f += 2;
    }

    face_data_f[mr_f] = 0.5 * (face_data_c[mr_c] + face_data_c[mr_c + i_rowsize_coarse]);
    face_data_f[mr_f-1] = 0.5 * (face_data_c[mr_c] + face_data_c[mr_c + i_rowsize_coarse - 1]);
    face_data_f[mr_f + rowsize_fine - 1 - 1] = 0.5 * (face_data_c[mr_c + i_rowsize_coarse] + face_data_c[mr_c + i_rowsize_coarse - 1]);

    mr_c += 3;
    mr_f += rowsize_fine - 1 + 3;

    rowsize_fine -= 2;
    i_rowsize_coarse -= 1;
  }
  double sum = 0;
  for(int i = 0; i < num_microvertices_per_face(level-2+1)*sizeof(double); ++i){
    sum += face_data_f[i];
  }
  std::cout << sum << std::endl;
}



int main()
{
  prolongate(6);
  return 0;
}