#pragma once
#include <array>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
namespace etch {
// Roots on a bounded segment. Derivative roots split each polynomial into monotone intervals.
inline double polynomial_value(const std::vector<double>& coefficients,double t) {
    double value=0;for(auto i=coefficients.rbegin();i!=coefficients.rend();++i)value=value*t+*i;return value;
}
inline std::vector<double> segment_roots(std::vector<double> coefficients) {
    double scale=0;for(double c:coefficients)scale+=std::abs(c);
    if(scale==0)return {};
    while(coefficients.size()>1&&std::abs(coefficients.back())<1e-14*scale)coefficients.pop_back();
    int degree=int(coefficients.size())-1;
    if(degree==0)return {};
    if(degree==1) {
        double t=-coefficients[0]/coefficients[1];
        // Adjacent DDA intervals share endpoints; retain roots displaced by rounding.
        return t>=-1e-13&&t<=1+1e-13?std::vector<double>{std::clamp(t,0.0,1.0)}:std::vector<double>{};
    }
    std::vector<double> derivative(degree);for(int i=1;i<=degree;++i)derivative[i-1]=i*coefficients[i];
    std::vector<double> intervals{0};
    for(double t:segment_roots(derivative))if(t>0&&t<1)intervals.push_back(t);
    intervals.push_back(1);std::vector<double> roots;
    for(double t:intervals)if(std::abs(polynomial_value(coefficients,t))<=1e-13*scale)roots.push_back(t);
    for(size_t i=1;i<intervals.size();++i) {
        double left=intervals[i-1],right=intervals[i],fl=polynomial_value(coefficients,left),fr=polynomial_value(coefficients,right);
        if(fl==0||fr==0||(fl<0)==(fr<0))continue;
        for(int iteration=0;iteration<55;++iteration) {
            double middle=(left+right)*0.5,fm=polynomial_value(coefficients,middle);
            if((fl<0)==(fm<0)){left=middle;fl=fm;}else right=middle;
        }
        roots.push_back((left+right)*0.5);
    }
    std::sort(roots.begin(),roots.end());
    roots.erase(std::unique(roots.begin(),roots.end(),[](double a,double b){return std::abs(a-b)<1e-12;}),roots.end());return roots;
}

// Compact partition of unity. Each jet contains distinct-axis mixed derivatives.
// The geometric window reuses the local fitting support; it is independent of cell-removal weights.
template<int D> struct ImplicitBlend {
    using Node=std::array<double,1<<D>;
    Node weights{},numerator{};
    void add(const std::array<double,D>& offset,const std::array<double,D>& radius,const Node& local) {
        std::array<std::array<double,2>,D> factors{};
        for(int axis=0;axis<D;++axis) {
            double t=offset[axis]/radius[axis];if(std::abs(t)>=1)return;
            double a=1-t*t;factors[axis]={a*a*a*a,-8*t*a*a*a/radius[axis]};
        }
        Node weight{};
        for(int mask=0;mask<(1<<D);++mask) {
            weight[mask]=1;for(int axis=0;axis<D;++axis)weight[mask]*=factors[axis][(mask>>axis)&1];
            weights[mask]+=weight[mask];
        }
        for(int mask=0;mask<(1<<D);++mask)for(int subset=mask;;subset=(subset-1)&mask) {
            numerator[mask]+=weight[subset]*local[mask^subset];if(subset==0)break;
        }
    }
    Node normalized() const {
        Node result{};result[0]=numerator[0]/weights[0];
        for(int mask=1;mask<(1<<D);++mask) {
            double value=numerator[mask];
            for(int subset=mask;subset;subset=(subset-1)&mask)value-=weights[subset]*result[mask^subset];
            result[mask]=value/weights[0];
        }
        return result;
    }
};

// Tensor Hermite reconstruction: node[mask] stores mixed derivatives in physical coordinates.
// Shared node values and derivatives make the field C1 across reconstruction-cell faces.
template<int D> struct ImplicitCell {
    static constexpr int size=1<<(2*D);
    using Node=std::array<double,1<<D>;
    std::array<double,size> coefficients{};
    double minimum=0,maximum=0;
    static std::array<double,size> transform(std::array<double,size> data,const double (&matrix)[4][4]) {
        int stride=1;
        for(int axis=0;axis<D;++axis,stride*=4) {
            auto previous=data;
            for(int i=0;i<size;++i) {
                int digit=(i/stride)%4,base=i-digit*stride;data[i]=0;
                for(int j=0;j<4;++j)data[i]+=matrix[digit][j]*previous[base+j*stride];
            }
        }
        return data;
    }
    ImplicitCell()=default;
    ImplicitCell(const std::array<Node,1<<D>& nodes,double spacing) {
        std::array<double,size> data{};
        for(int corner=0;corner<(1<<D);++corner)for(int mask=0;mask<(1<<D);++mask) {
            int index=0,stride=1;double scale=1;
            for(int axis=0;axis<D;++axis,stride*=4) {
                index+=(((corner>>axis)&1)+2*((mask>>axis)&1))*stride;
                if(mask&(1<<axis))scale*=spacing;
            }
            data[index]=nodes[corner][mask]*scale;
        }
        const double hermite[4][4]={{1,0,0,0},{0,0,1,0},{-3,3,-2,-1},{2,-2,1,1}};
        coefficients=transform(data,hermite);
        const double bernstein[4][4]={{1,0,0,0},{1,1.0/3,0,0},{1,2.0/3,1.0/3,0},{1,1,1,1}};
        auto bounds=transform(coefficients,bernstein);
        auto range=std::minmax_element(bounds.begin(),bounds.end());
        double scale=0;for(double c:coefficients)scale+=std::abs(c);
        double rounding=64*D*std::numeric_limits<double>::epsilon()*scale;
        minimum=*range.first-rounding;maximum=*range.second+rounding;
    }
    double value(const std::array<double,D>& point) const {
        std::array<std::array<double,4>,D> powers{};
        for(int axis=0;axis<D;++axis){powers[axis][0]=1;for(int k=1;k<4;++k)powers[axis][k]=powers[axis][k-1]*point[axis];}
        double result=0;
        for(int i=0;i<size;++i){double term=coefficients[i];for(int axis=0,index=i;axis<D;++axis,index/=4)term*=powers[axis][index%4];result+=term;}
        return result;
    }
    std::array<double,D> gradient(const std::array<double,D>& point) const {
        std::array<std::array<double,4>,D> powers{};std::array<double,D> result{};
        for(int axis=0;axis<D;++axis){powers[axis][0]=1;for(int k=1;k<4;++k)powers[axis][k]=powers[axis][k-1]*point[axis];}
        for(int i=0;i<size;++i)for(int derivative=0;derivative<D;++derivative) {
            double term=coefficients[i];
            for(int axis=0,index=i;axis<D;++axis,index/=4) {
                int degree=index%4;if(axis==derivative){if(degree==0){term=0;break;}term*=degree*powers[axis][degree-1];}else term*=powers[axis][degree];
            }
            result[derivative]+=term;
        }
        return result;
    }
    std::pair<double,double> last_axis_range(double coordinate) const {
        constexpr int count=size/4;std::array<double,count> data{};
        for(int i=0;i<count;++i)for(int k=3;k>=0;--k)data[i]=data[i]*coordinate+coefficients[i+k*count];
        const double bernstein[4][4]={{1,0,0,0},{1,1.0/3,0,0},{1,2.0/3,1.0/3,0},{1,1,1,1}};
        int stride=1;
        for(int axis=0;axis<D-1;++axis,stride*=4) {
            auto previous=data;
            for(int i=0;i<count;++i) {
                int digit=(i/stride)%4,base=i-digit*stride;data[i]=0;
                for(int j=0;j<4;++j)data[i]+=bernstein[digit][j]*previous[base+j*stride];
            }
        }
        auto bounds=std::minmax_element(data.begin(),data.end());
        double scale=0;for(double c:coefficients)scale+=std::abs(c);
        double rounding=64*D*std::numeric_limits<double>::epsilon()*scale;
        return {*bounds.first-rounding,*bounds.second+rounding};
    }
    std::vector<double> intersections(const std::array<double,D>& start,const std::array<double,D>& end) const {
        if(minimum>0||maximum<0)return {};
        std::array<std::array<std::array<double,4>,4>,D> powers{};
        for(int axis=0;axis<D;++axis) {
            powers[axis][0][0]=1;double delta=end[axis]-start[axis];
            for(int degree=1;degree<4;++degree)for(int k=0;k<=degree;++k)
                powers[axis][degree][k]=powers[axis][degree-1][k]*start[axis]+(k?powers[axis][degree-1][k-1]*delta:0);
        }
        std::vector<double> polynomial(3*D+1);
        for(int i=0;i<size;++i) {
            std::array<double,3*D+1> product{};product[0]=coefficients[i];int degree=0;
            for(int axis=0,index=i;axis<D;++axis,index/=4) {
                int power=index%4;auto previous=product;product.fill(0);
                for(int k=0;k<=degree;++k)for(int j=0;j<=power;++j)product[k+j]+=previous[k]*powers[axis][power][j];
                degree+=power;
            }
            for(int k=0;k<=degree;++k)polynomial[k]+=product[k];
        }
        return segment_roots(polynomial);
    }
};
}
