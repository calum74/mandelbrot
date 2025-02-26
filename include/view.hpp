namespace fractals
{
    // Manages the multithreaded animation and rendering
    class view
    {
    public:
        pixmap<error_value<double>> values;

        void set_size(int w, int h);
        void set_calculation(std::shared_ptr<fractal_calculation> & calculation);

        class listener
        {
        public:
            virtual void redraw() =0;
            virtual void update_status(const rendering_metrics&);
        };
    
    };
}