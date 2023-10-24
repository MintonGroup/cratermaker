struct surface_type {
        double *elev_data;
        int elev_shape[2];
};
double c_double;
extern struct surface_type* bind_surface_init(int gridsize);
extern void bind_surface_final(struct surface_type* obj);