from beluga.visualization import BelugaPlot

plots = BelugaPlot('./data.dill',default_sol=-1,default_step=-1)

plots.add_plot().line('v','h')                    \
                .xlabel('v(t)').ylabel('h(t)')      \
                .title('Trajectory')
'''
plots.add_plot().line('t','theta')                    \
                .xlabel('t (s)').ylabel('theta (degrees)')      \
                .title('Control history')

plots.add_plot().line_series('x','y') \
                .xlabel('x(t)').ylabel('y(t)')      \
                .title('Trajectory history')
'''
plots.render()

#/Users/hsaranat/anaconda/lib/python3.5/site-packages/numexpr/necompiler.py
