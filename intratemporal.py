#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 16:50:59 2020

@author: egorkozlov
"""
import numpy as np
from scipy.optimize import brentq

def int_with_x(m_in,newton=True,step=1e-6,nint=2000,*,A,alp,sig,xi,lam,kap,lbr,qlb=0.0,cleft=0.05):
    m_in = np.atleast_1d(m_in)    
    
    def foc_expression(x):
        return ((A/alp)**(1/sig))*(x**lam + kap*(1-lbr)**lam)**((lam+xi-1)/(lam*sig)) * \
                x**((1-lam)/sig) + x
    
    
    
    rf_res = lambda x : foc_expression(x) - m_in[0]
    
    x0 = 1e-6
    
    factor = 1.0
    
    for ib in range(15):
        try:
            x0 = brentq(rf_res,step,m_in[0]-step)
            break
        except:
            if ib == 14: print('brent failed...')
            step = step/10
            factor = 0.2 # expect slow convergence
            pass
        
    
    
    xgrid = np.linspace(0.5*x0,m_in.max(),nint)
    
    
    m_implied = foc_expression(xgrid)
          
    x_interpolated = np.maximum(np.interp(m_in,m_implied,xgrid),x0)
    
    xout = x_interpolated
    if newton: 
        # perform  iteration of Newton method
        def foc_deriv(x):
            logder = ((lam+xi-1)/(sig)) * (x**(lam-1))/(x**lam + kap*(1-lbr)**lam) + \
                ((1-lam)/sig)/x
            return 1 + (foc_expression(x) - x)*logder
        
        
        for i in range(10):
            x_start = xout
            f_res_start = foc_expression(x_start) - m_in
            f_der = foc_deriv(x_start)
    
            xout = np.maximum( x_start - factor*f_res_start/f_der, x0 )
    
            f_res = foc_expression(xout) - m_in
            
            assert not np.any(np.isnan(f_res))
            
            #print('At i = {} max residuals are {}'.format(i,np.abs(f_res).max()))
            
            if np.abs(f_res).max() < 1e-4:
                break
            if np.abs(f_res).max() > np.abs(f_res_start).max() and i > 2:
                print('Intratemporal: backing up...')
                imax = np.abs(f_res).argmax()
                print('Max residual is {} at x = {}'.
                      format(f_res[imax],x_start[imax]))                
                i_worse = (np.abs(f_res) > np.abs(f_res_start))
                
                xout[i_worse] = x_start[i_worse]
                break
        
        
    
    
    
    
    def qfun(x,l):
        return (x**lam + kap*(1-l)**lam)**(1/lam)
    
    def x_inv(q,l):
        return (q**lam - kap*(1-l)**lam)**(1/lam)
    
    def util(c,q):
        uc = A*c**(1-sig)/(1-sig)
        ux = alp*q**(1-xi)/(1-xi)
        return uc + ux
    
    
    
    cout = m_in - xout
    qout = qfun(xout,lbr)
    uout = util(cout,qout)
    
    
    
    # enforcing lower bound
    
    if qlb > 0:
        # enforcing the bounds
        q_min = qfun((1-cleft)*m_in,lbr)
        qlb_fixed = np.minimum( q_min, qlb )
        to_fix = (qout < qlb)
        xout[to_fix] = x_inv(qlb_fixed[to_fix],lbr)
        qout[to_fix] = qfun(xout[to_fix],lbr)
        cout[to_fix] = m_in[to_fix] - xout[to_fix]
        uout[to_fix] = util(cout[to_fix],qout[to_fix])
        assert np.all(qout >= qlb_fixed-1e-5)
    
    return xout, cout, uout, qout



def int_no_x(m_in,nint=2000,*,A,sig):
    cout = m_in
    uout = A*(cout**(1-sig))/(1-sig)
    xout = np.zeros(cout.shape,dtype=cout.dtype)
    qout = np.zeros(cout.shape,dtype=cout.dtype)
    
    return xout, cout, uout, qout
    
    