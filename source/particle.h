/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Base class for particle systems
 *
 ******************************************************************************/

#ifndef _PARTICLE_H
#define _PARTICLE_H

#include <vector>
#include "grid.h"
#include "vectorbase.h"
#include "integrator.h"
#include "randomstream.h"
namespace Manta {

// fwd decl
template<class T> class Grid;
class ParticleDataBase;
template<class T> class ParticleDataImpl;

//! Baseclass for particle systems. Does not implement any data
PYTHON() class ParticleBase : public PbClass {
public:
	enum SystemType { BASE=0, PARTICLE, VORTEX, FILAMENT, FLIP, TURBULENCE, INDEX };
	
	enum ParticleStatus {
		PNONE         = 0,
		PNEW          = (1<<0),  // particles newly created in this step
		PSPRAY        = (1<<1),  // secondary particle types
		PBUBBLE       = (1<<2),
		PFOAM         = (1<<3),
		PTRACER       = (1<<4),
		PDELETE       = (1<<10), // mark as deleted, will be deleted in next compress() step
		PINVALID      = (1<<30), // unused
	};

	PYTHON() ParticleBase(FluidSolver* parent);
	virtual ~ParticleBase();

	//! copy all the particle data thats registered with the other particle system to this one
	virtual void cloneParticleData(ParticleBase* nm);

	virtual SystemType getType() const { return BASE; }
	virtual std::string infoString() const; 
	virtual ParticleBase* clone() { assertMsg( false , "Dont use, override..."); return NULL; } 

	//! slow virtual function to query size, do not use in kernels! use size() instead
	virtual IndexInt getSizeSlow() const { assertMsg( false , "Dont use, override..."); return 0; }

	//! add a position as potential candidate for new particle (todo, make usable from parallel threads)
	// inline void addBuffered(const Vec3& pos, int flag=0);
    inline void addBuffered(const Vec3& pos, int flag=0) {
	    mNewBufferPos.push_back(pos);
	    mNewBufferFlag.push_back(flag);
    }

	//! particle data functions

	//! create a particle data object
	PYTHON() PbClass* create(PbType type, PbTypeVec T=PbTypeVec(), const std::string& name = "");
	//! add a particle data field, set its parent particle-system pointer
	void registerPdata(ParticleDataBase* pdata);
	void registerPdataReal(ParticleDataImpl<Real>* pdata);
	void registerPdataVec3(ParticleDataImpl<Vec3>* pdata);
	void registerPdataInt (ParticleDataImpl<int >* pdata);
	//! remove a particle data entry
	void deregister(ParticleDataBase* pdata);
	//! add one zero entry to all data fields
	void addAllPdata();
	// note - deletion of pdata is handled in compress function

	//! how many are there?
	IndexInt getNumPdata() const { return mPartData.size(); }
	//! access one of the fields
	ParticleDataBase* getPdata(int i) { return mPartData[i]; }

protected:  
	//! new particle candidates
	std::vector<Vec3> mNewBufferPos;
	std::vector<int> mNewBufferFlag;

	//! allow automatic compression / resize? disallowed for, eg, flip particle systems
	bool mAllowCompress;

	//! store particle data , each pointer has its own storage vector of a certain type (int, real, vec3)
	std::vector<ParticleDataBase*> mPartData;
	//! lists of different types, for fast operations w/o virtual function calls (all calls necessary per particle)
	std::vector< ParticleDataImpl<Real> *> mPdataReal;
	std::vector< ParticleDataImpl<Vec3> *> mPdataVec3;
	std::vector< ParticleDataImpl<int> *>  mPdataInt;
	//! indicate that pdata of this particle system is copied, and needs to be freed
	bool mFreePdata;
};


//! Main class for particle systems
/*! Basetype S must at least contain flag, pos fields */
PYTHON() template<class S> class ParticleSystem : public ParticleBase {
public:    
	PYTHON() ParticleSystem(FluidSolver* parent) : ParticleBase(parent), mDeletes(0), mDeleteChunk(0) {}
	virtual ~ParticleSystem() {};
	
	virtual SystemType getType() const { return S::getType(); };
	
	//! accessors
	inline S& operator[](IndexInt idx)             { DEBUG_ONLY(checkPartIndex(idx)); return mData[idx]; }
	inline const S& operator[](IndexInt idx) const { DEBUG_ONLY(checkPartIndex(idx)); return mData[idx]; }
	//! return size of container
	//! note , python binding disabled for now! cannot yet deal with long-long types
	inline IndexInt size() const { return mData.size(); }
	//! slow virtual function of base class, also returns size
	virtual IndexInt getSizeSlow() const { return size(); }
	//! note , special call for python, note - doesnt support more than 2b parts!
	PYTHON() int pySize() const { return (int)mData.size(); }

	//! query status
	inline int  getStatus(IndexInt idx) const { DEBUG_ONLY(checkPartIndex(idx)); return mData[idx].flag; }
	inline bool isActive(IndexInt idx) const { DEBUG_ONLY(checkPartIndex(idx)); return (mData[idx].flag & PDELETE) == 0; }
	inline bool isSpray(IndexInt idx) const { DEBUG_ONLY(checkPartIndex(idx)); return (mData[idx].flag & PSPRAY); }
	inline bool isBubble(IndexInt idx) const { DEBUG_ONLY(checkPartIndex(idx)); return (mData[idx].flag & PBUBBLE); }
	inline bool isFoam(IndexInt idx) const { DEBUG_ONLY(checkPartIndex(idx)); return (mData[idx].flag & PFOAM); }
	inline bool isTracer(IndexInt idx) const { DEBUG_ONLY(checkPartIndex(idx)); return (mData[idx].flag & PTRACER); }

	//! update status
	inline void setStatus(IndexInt idx, const int status) { DEBUG_ONLY(checkPartIndex(idx)); mData[idx].flag = status; }
	
	//! safe accessor for python
	PYTHON() void setPos(const IndexInt idx, const Vec3& pos) { DEBUG_ONLY(checkPartIndex(idx)); mData[idx].pos = pos; }
	PYTHON() const Vec3& getPos(const IndexInt idx) const     { DEBUG_ONLY(checkPartIndex(idx)); return mData[idx].pos; }
	//! copy all positions into pdata vec3 field
	PYTHON() void getPosPdata(ParticleDataImpl<Vec3>& target) const;
	PYTHON() void setPosPdata(const ParticleDataImpl<Vec3>& source);
	//! transform coordinate system from one grid size to another (usually upon load)
	void transformPositions( Vec3i dimOld, Vec3i dimNew );

	//! explicitly trigger compression from outside
	void doCompress() { if ( mDeletes > mDeleteChunk) compress(); }
	//! insert buffered positions as new particles, update additional particle data
	void insertBufferedParticles();
	//! resize data vector, and all pdata fields
	void resizeAll(IndexInt newsize);
	
	//! adding and deleting
	inline void kill(IndexInt idx);
	IndexInt add(const S& data);
	//! remove all particles, init 0 length arrays (also pdata)
	PYTHON() void clear();
			
	//! Advect particle in grid velocity field
	PYTHON() void advectInGrid(const FlagGrid &flags, const MACGrid &vel, const int integrationMode, const bool deleteInObstacle=true, const bool stopInObstacle=true, const bool skipNew=false, const ParticleDataImpl<int> *ptype=NULL, const int exclude=0);
	
	//! Project particles outside obstacles
	PYTHON() void projectOutside(Grid<Vec3> &gradient);
	PYTHON() void projectOutOfBnd(const FlagGrid &flags, const Real bnd, const std::string &plane="xXyYzZ", const ParticleDataImpl<int> *ptype=NULL, const int exclude=0);
	
	virtual ParticleBase* clone();
	virtual std::string infoString() const;

	//! debugging
	inline void checkPartIndex(IndexInt idx) const;
	
protected:  
	//! deletion count , and interval for re-compressing 
	IndexInt mDeletes, mDeleteChunk;
	//! the particle data
	std::vector<S> mData;    

	//! reduce storage , called by doCompress
	virtual void compress(); 
};

//******************************************************************************

//! Simplest data class for particle systems
//! contains a position and an int flag; note that these are deprectated, and will at
//! some point be replaced by the more flexible pdata fields. For now manually copy with
//! getPosPdata / setPosPdata.
struct BasicParticleData {
public:
	BasicParticleData() : pos(0.), flag(0) {}
	BasicParticleData(const Vec3& p) : pos(p), flag(0) {}
	static ParticleBase::SystemType getType() { return ParticleBase::PARTICLE; }

	//! data (note, this size is currently hard coded for uni i/o)
	Vec3 pos;
	int  flag;
};

PYTHON() class BasicParticleSystem : public ParticleSystem<BasicParticleData> {
public:
	PYTHON() BasicParticleSystem(FluidSolver* parent);
	
	//! file io
	PYTHON() int save(const std::string name);
	PYTHON() int load(const std::string name);

	//! save to text file
	void writeParticlesText(const std::string name) const;
	//! other output formats
	void writeParticlesRawPositionsGz(const std::string name) const;
	void writeParticlesRawVelocityGz(const std::string name) const;

	//! read from other particle system (with resize)
	PYTHON() void readParticles(BasicParticleSystem* from);

	//! add particles in python
	PYTHON() void addParticle(Vec3 pos) { add(BasicParticleData(pos)); }

	//! dangerous, get low level access - avoid usage, only used in vortex filament advection for now
	std::vector<BasicParticleData>& getData() { return mData; }

	PYTHON() void printParts(IndexInt start=-1, IndexInt stop=-1, bool printIndex=false);

	//! get data pointer of particle data
	PYTHON() std::string getDataPointer();
};


//******************************************************************************

//! Index into other particle system
//  used for grid based neighborhood searches on generic particle systems (stores
//  only active particles, and reduces copied data)
//  note - pos & flag are disabled here, do not use!
struct ParticleIndexData {
public:
	ParticleIndexData() : sourceIndex(0) {}
	static ParticleBase::SystemType getType() { return ParticleBase::INDEX; }

	IndexInt  sourceIndex; // index of this particle in the original particle system
	//! note - the following two are needed for template instantiation, but not used
	//! for the particle index system (use values from original one!)
	static Vec3 pos;  // do not use... 
	static int  flag; // not needed usally 
	//Vec3 pos; // enable for debugging
};

PYTHON() class ParticleIndexSystem : public ParticleSystem<ParticleIndexData> {
public:
	PYTHON() ParticleIndexSystem(FluidSolver* parent) : ParticleSystem<ParticleIndexData>(parent) {};
	
	//! we only need a resize function...
	void resize(IndexInt size) { mData.resize(size); }
};



//******************************************************************************

//! Particle set with connectivity
PYTHON() template<class DATA, class CON>
class ConnectedParticleSystem : public ParticleSystem<DATA> {
public:
	PYTHON() ConnectedParticleSystem(FluidSolver* parent) : ParticleSystem<DATA>(parent) {}
	
	//! accessors
	inline bool isSegActive(int i) { return (mSegments[i].flag & ParticleBase::PDELETE) == 0; }    
	inline int segSize() const { return mSegments.size(); }    
	inline CON& seg(int i) { return mSegments[i]; }
	inline const CON& seg(int i) const { return mSegments[i]; }
		
	virtual ParticleBase* clone();
	
protected:
	std::vector<CON> mSegments;
	virtual void compress();    
};

//******************************************************************************

//! abstract interface for particle data
PYTHON() class ParticleDataBase : public PbClass {
public:
	PYTHON() ParticleDataBase(FluidSolver* parent);
	virtual ~ParticleDataBase(); 

	//! data type IDs, in line with those for grids
	enum PdataType { TypeNone = 0, TypeReal = 1, TypeInt = 2, TypeVec3 = 4 };

	//! interface functions, using assert instead of pure virtual for python compatibility
	virtual IndexInt  getSizeSlow() const { assertMsg( false , "Dont use, override..."); return 0; } 
	virtual void addEntry()   { assertMsg( false , "Dont use, override..."); return;   }
	virtual ParticleDataBase* clone() { assertMsg( false , "Dont use, override..."); return NULL; }
	virtual PdataType getType() const { assertMsg( false , "Dont use, override..."); return TypeNone; } 
	virtual void resize(IndexInt size)     { assertMsg( false , "Dont use, override..."); return;  }
	virtual void copyValueSlow(IndexInt from, IndexInt to) { assertMsg( false , "Dont use, override..."); return;  }

	//! set / get base pointer to parent particle system
	void setParticleSys(ParticleBase* set) { mpParticleSys = set; }
	ParticleBase* getParticleSys()         { return mpParticleSys; }

	//! debugging
	inline void checkPartIndex(IndexInt idx) const;

protected:
	ParticleBase* mpParticleSys;
};


//! abstract interface for particle data
PYTHON() template<class T>
class ParticleDataImpl : public ParticleDataBase {
public:
	PYTHON() ParticleDataImpl(FluidSolver* parent);
	ParticleDataImpl(FluidSolver* parent, ParticleDataImpl<T>* other);
	virtual ~ParticleDataImpl();

	//! access data
	inline       T& get(const IndexInt idx)              { DEBUG_ONLY(checkPartIndex(idx)); return mData[idx]; }
	inline const T& get(const IndexInt idx) const        { DEBUG_ONLY(checkPartIndex(idx)); return mData[idx]; }
	inline       T& operator[](const IndexInt idx)       { DEBUG_ONLY(checkPartIndex(idx)); return mData[idx]; }
	inline const T& operator[](const IndexInt idx) const { DEBUG_ONLY(checkPartIndex(idx)); return mData[idx]; }

	//! set data
	inline       void set(const IndexInt idx, T& val)    { DEBUG_ONLY(checkPartIndex(idx)); mData[idx] = val; }

	//! set all values to 0, note - different from particleSystem::clear! doesnt modify size of array (has to stay in sync with parent system)
	PYTHON() void clear();

	//! set grid from which to get data...
	PYTHON() void setSource(Grid<T>* grid, bool isMAC=false );

	//! particle data base interface
	virtual IndexInt  getSizeSlow() const;
	virtual void addEntry();
	virtual ParticleDataBase* clone();
	virtual PdataType getType() const;
	virtual void resize(IndexInt s);
	virtual void copyValueSlow(IndexInt from, IndexInt to);

	IndexInt  size() const { return mData.size(); }

	//! fast inlined functions for per particle operations
	inline void copyValue(IndexInt from, IndexInt to) { get(to) = get(from); }
	void initNewValue(IndexInt idx, Vec3 pos);

	//! python interface (similar to grid data)
	PYTHON() ParticleDataImpl<T>& copyFrom(const ParticleDataImpl<T>& a);
	PYTHON() void setConst(const T &s);
	PYTHON() void setConstRange(const T &s, const int begin, const int end);
	PYTHON() void add(const ParticleDataImpl<T>& a);
	PYTHON() void sub(const ParticleDataImpl<T>& a);
	PYTHON() void addConst(const T &s);
	PYTHON() void addScaled(const ParticleDataImpl<T>& a, const T& factor);
	PYTHON() void mult(const ParticleDataImpl<T>& a);
	PYTHON() void multConst(const T &s);
	PYTHON() void safeDiv(const ParticleDataImpl<T>& a);
	PYTHON() void clamp(const Real vmin, const Real vmax);
	PYTHON() void clampMin(const Real vmin);
	PYTHON() void clampMax(const Real vmax);

	PYTHON() Real getMaxAbs() const;
	PYTHON() Real getMax() const;
	PYTHON() Real getMin() const;

	PYTHON() T    sum(const ParticleDataImpl<int> *t=NULL, const int itype=0) const;
	PYTHON() Real sumSquare() const;
	PYTHON() Real sumMagnitude() const;

	//! special, set if int flag in t has "flag"
	PYTHON() void setConstIntFlag(const T &s, const ParticleDataImpl<int> &t, const int flag);

	PYTHON() void printPdata(IndexInt start=-1, IndexInt stop=-1, bool printIndex=false);
	
	//! file io
	PYTHON() int save(const std::string name);
	PYTHON() int load(const std::string name);

	//! get data pointer of particle data
	PYTHON() std::string getDataPointer();
protected:
	//! data storage
	std::vector<T> mData; 

	//! optionally , we might have an associated grid from which to grab new data
	Grid<T>* mpGridSource;
	//! unfortunately , we need to distinguish mac vs regular vec3
	bool mGridSourceMAC;
};

PYTHON() alias ParticleDataImpl<int>  PdataInt;
PYTHON() alias ParticleDataImpl<Real> PdataReal;
PYTHON() alias ParticleDataImpl<Vec3> PdataVec3;




//! count by type flag
int countParticles(const ParticleDataImpl<int> &t, const int flag);
// set contents to zero, as for a grid
template<class T>
void ParticleDataImpl<T>::clear() {
	for(IndexInt i=0; i<(IndexInt)mData.size(); ++i) mData[i] = 0.;
}
} // namespace

#endif

