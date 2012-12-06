#ifndef MERGEDESCRIPTOR_HPP_
#define MERGEDESCRIPTOR_HPP_

#include "pctg/BestCtgAlignment.hpp"
#include "assembly/contig.hpp"

typedef enum { LINEAR_MERGE, FORK_MERGE } __attribute__((packed)) MergeType;
typedef enum { MIS_MASTER, MIS_SLAVE, REPEAT, UNKNOWN } __attribute__((packed)) ForkType;

struct MergeBlock
{
	uint64_t vertex;

	int32_t m_id;
	int32_t m_start;
	int32_t m_end;

	int32_t s_id;
	int32_t s_start;
	int32_t s_end;
    
    bool valid;

	bool align_rev;
	bool align_ok;

	// relative to orginal strand of contigs
	bool m_ltail;
	bool m_rtail;
	bool s_ltail;
	bool s_rtail;

	// relative to order of merge
	bool ext_slave_next;
	bool ext_slave_prev;

	bool m_rev;
	bool s_rev;
};

struct MergeDescriptor
{
    MergeType mergeType;
    ForkType forkType;

	int32_t masterId;
	int32_t slaveId;

	std::pair<int64_t,int64_t> masterMainAlign;
	std::pair<int64_t,int64_t> slaveMainAlign;

    std::pair<int64_t,int64_t> masterPos;
	std::pair<int64_t,int64_t> slavePos;

	std::pair<int64_t,int64_t> masterBlocks;
	std::pair<int64_t,int64_t> slaveBlocks;

	BestCtgAlignment *align;
};

typedef std::list< std::list<MergeDescriptor> > MergeDescriptorLists;
typedef std::list< std::list<MergeBlock> > MergeBlockLists;


struct MergeStruct
{
	int32_t id;
	Contig *ctg;
	bool reversed;
	int64_t pos;

	std::pair<int64_t,int64_t> limits;
};

#endif //MERGEDESCRIPTOR_HPP_
