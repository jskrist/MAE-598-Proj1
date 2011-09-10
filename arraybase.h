/*********************************************
Utility Library Function
Copyright(c) 2006, S. D. Rajan
All rights reserved

Object-Oriented Numerical Analysis
*********************************************/
#ifndef __RAJAN_ARRAYBASE_H__
#define __RAJAN_ARRAYBASE_H__

class CArrayBase
{
    public:
        CArrayBase ();
        ~CArrayBase ();

        void static ShowStatistics ();

    protected:
        static double m_dAllocated;
        static double m_dDeAllocated;
};

#endif

