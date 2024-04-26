
typedef struct c_PerlinArguments {
      double damp;
      double damp0;
      double damp_scale;
      double freq;
      double gain;
      double gain0;
      double lacunarity;
      double noise_height;
      double pers;
      double slope;
      double warp;
      double warp0;
}PerlinArguments;

double bind_realistic_perlin_noise_one(const char *model, double x, double y, double z, int num_octaves, double *anchor, double damp, double damp0, double damp_scale, double freq, double gain, double gain0, double lacunarity, double noise_height, double pers, double slope, double warp, double warp0);

void bind_realistic_perlin_noise_all(const char *model, double *x, double *y, double *z, int num_elements, int num_octaves, double *anchor, double damp, double damp0, double damp_scale, double freq, double gain, double gain0, double lacunarity, double noise_height, double pers, double slope, double warp, double warp0, double *noise);