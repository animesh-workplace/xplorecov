import numpy as np
from .models import LLMCache
from django.db.models import F
from django.conf import settings
from django.utils import timezone
from typing import Optional, Tuple
from asgiref.sync import async_to_sync
from scipy.spatial.distance import cosine
import json, hashlib, pickle, redis, asyncio
from sentence_transformers import SentenceTransformer


class LLMCacheManager:
    def __init__(self):
        self.redis_client = redis.Redis(
            db=0,
            decode_responses=True,
            host="127.0.0.1",
            port="6379",
        )
        self.cache_timeout = 3600  # 1 hour Redis cache timeout
        # Initialize the sentence transformer model
        self.sentence_transformer = SentenceTransformer("all-MiniLM-L6-v2")
        self.similarity_threshold = 0.85  # Adjust this threshold based on your needs

    def _generate_cache_key(self, query: str) -> str:
        """
        Generate a unique and consistent cache key for a query.

        Args:
            query (str): The user's query string

        Returns:
            str: A cache key string in the format 'llm_cache:{md5_hash}'
        """
        # Normalize the query by lowercasing and removing extra whitespace
        normalized_query = " ".join(query.lower().split())

        # Generate MD5 hash of the normalized query
        query_hash = hashlib.md5(normalized_query.encode("utf-8")).hexdigest()

        # Create the cache key with a prefix for namespace isolation
        cache_key = f"llm_cache:{query_hash}"

        return cache_key

    def _store_in_redis(
        self, cache_key: str, cache_data: dict, timeout: Optional[int] = None
    ) -> None:
        """
        Store data in Redis with optional custom timeout.

        Args:
            cache_key (str): The key to store the data under
            cache_data (dict): The data to store
            timeout (Optional[int]): Custom timeout in seconds, defaults to self.cache_timeout
        """
        if timeout is None:
            timeout = self.cache_timeout

        try:
            serialized_data = json.dumps(cache_data)
            self.redis_client.setex(cache_key, timeout, serialized_data)
        except (redis.RedisError, json.JSONDecodeError) as e:
            # Log the error but don't raise it to prevent breaking the application
            print(f"Redis storage error: {str(e)}")

    async def _update_hit_count_async(self, user_query: str) -> None:
        """
        Asynchronously update the hit count for a cache entry.

        Args:
            user_query (str): The original user query to update
        """
        try:
            # Use F() to prevent race conditions
            await LLMCache.objects.filter(user_query=user_query).aupdate(
                hit_count=F("hit_count") + 1, last_accessed=timezone.now()
            )
        except Exception as e:
            # Log the error but don't raise it
            print(f"Hit count update error: {str(e)}")

    def _update_hit_count(self, user_query: str) -> None:
        """
        Synchronous wrapper for hit count update.

        Args:
            user_query (str): The original user query to update
        """
        try:
            # For Django sync views, use sync_to_async
            async_to_sync(self._update_hit_count_async)(user_query)
        except Exception as e:
            # Fallback to synchronous update if async fails
            try:
                LLMCache.objects.filter(user_query=user_query).update(
                    hit_count=F("hit_count") + 1, last_accessed=timezone.now()
                )
            except Exception as fallback_e:
                print(f"Hit count update error (fallback): {str(fallback_e)}")

    def _generate_embedding(self, query: str) -> np.ndarray:
        """Generate embedding for a query."""
        return self.sentence_transformer.encode([query], normalize_embeddings=True)[0]

    def _calculate_similarity(
        self, embedding1: np.ndarray, embedding2: np.ndarray
    ) -> float:
        """Calculate cosine similarity between two embeddings."""
        print(embedding1, embedding2)
        return 1 - cosine(embedding1, embedding2)

    def find_similar_query(self, user_query: str) -> Optional[Tuple[LLMCache, float]]:
        """Find the most similar cached query above the similarity threshold."""
        query_embedding = self._generate_embedding(user_query)

        # Get all cached queries
        cached_queries = LLMCache.objects.all()
        max_similarity = 0
        most_similar_cache = None

        for cache_entry in cached_queries:
            cached_embedding = np.frombuffer(cache_entry.query_embedding)
            similarity = self._calculate_similarity(query_embedding, cached_embedding)

            if similarity > max_similarity and similarity >= self.similarity_threshold:
                max_similarity = similarity
                most_similar_cache = cache_entry

        if most_similar_cache:
            return most_similar_cache, max_similarity
        return None

    def get_cached_response(self, user_query: str) -> Optional[dict]:
        """Try to get cached response, including similar queries."""
        # First try exact match in Redis
        cache_key = self._generate_cache_key(user_query)
        redis_cache = self.redis_client.get(cache_key)

        if redis_cache:
            cache_data = json.loads(redis_cache)
            self._update_hit_count(user_query)
            return {**cache_data, "match_type": "exact"}

        # Then try exact match in SQLite
        try:
            cache_entry = LLMCache.objects.get(user_query=user_query)
            cache_data = self._prepare_cache_data(cache_entry)
            self._store_in_redis(cache_key, cache_data)
            return {**cache_data, "match_type": "exact"}
        except LLMCache.DoesNotExist:
            pass

        # Finally, try similarity matching
        similar_result = self.find_similar_query(user_query)
        if similar_result:
            cache_entry, similarity = similar_result
            cache_data = self._prepare_cache_data(cache_entry)
            return {
                **cache_data,
                "match_type": "similar",
                "similarity_score": similarity,
                "original_query": cache_entry.user_query,
            }

        return None

    def store_response(self, user_query: str, modified_query: str, generated_code: str):
        """Store response with query embedding."""
        cache_key = self._generate_cache_key(user_query)
        cache_data = {
            "modified_query": modified_query,
            "generated_code": generated_code,
        }

        # Generate and store embedding
        query_embedding = self._generate_embedding(user_query)

        # Store in Redis
        self.redis_client.setex(cache_key, self.cache_timeout, json.dumps(cache_data))

        # Store in SQLite with embedding
        LLMCache.objects.create(
            user_query=user_query,
            modified_query=modified_query,
            generated_code=generated_code,
            query_embedding=query_embedding.tobytes(),
        )

    def _prepare_cache_data(self, cache_entry: LLMCache) -> dict:
        """Prepare cache data from a cache entry."""
        return {
            "modified_query": cache_entry.modified_query,
            "generated_code": cache_entry.generated_code,
        }
