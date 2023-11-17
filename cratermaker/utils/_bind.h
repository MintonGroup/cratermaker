#define STRMAX 512

typedef struct c_target_body_type {
    double *elevation;
    int elevation_shape[2];
    char name[STRMAX];
}target_body_type;
extern struct c_target_body_type* bind_body_init(int ny, int nx);
extern void bind_body_final(struct c_target_body_type *obj);
extern void bind_body_set_name(struct c_target_body_type *obj, const char *c_string);
extern char* bind_body_get_name(struct c_target_body_type *obj);

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

double bind_perlin_noise(const char *model, double x, double y, double z, int num_octaves, double *anchor, double damp, double damp0, double damp_scale, double freq, double gain, double gain0, double lacunarity, double noise_height, double pers, double slope, double warp, double warp0);